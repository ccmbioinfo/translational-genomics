import argparse
from datetime import date
from json import JSONDecodeError
import pandas as pd
import requests

def get_pedigree_info(ped) -> [dict, str]:
    """
    JSON response contains pedigree and family objects
    The family object contains only those individuals in the pedigree that have been assigned Phenotips IDs
    The pedigree object contains members, an array of pedigree nodes where each node represents a family member,
    and the relationship object, an array describing the relationships between pedigree nodes (parent:child relationships)
    """
    members = (
        {}
    )  # create dictionary where keys are C4R IDs, and values are Phenotips ID, node ID, affected status, sex, and parental node IDs
    node_to_C4R = {}  # map pedigree node to C4R ID
    # iterate through familyMembers dict, rather than members dict, to retrieve only IDs of those with Phenotips ID (others are placeholders in pedigree)
    for member in ped["family"]["familyMembers"]:
        C4R = member["identifier"]
        pid = member["id"]  # Phenotips ID
        for member in ped["pedigree"]["members"]:
            properties = member.get("properties", None)
            if properties:
                if properties.get("id", None) == pid:
                    node_id = member["id"] # pedigree node id
                    sex = member["properties"].get("sex", None)
                    break
        if not sex:
            sex = get_sex(pid)

        node_to_C4R[node_id] = C4R
        # refer to members dict to retrieve affected status, if it exists
        affected = get_affected(pid, ped["pedigree"]["members"])
        # refer to relationships dict to retrieve parent node IDs
        parents = get_parents(node_id, ped["pedigree"]["relationships"])
        members[C4R] = {
            "pid": pid,
            "node_id": node_id,
            "affected": affected,
            "sex": sex,
            "parents": parents,
        }
    # get C4R IDs for parents
    for member in members:
        parents = members[member]["parents"]
        if parents:
            parents_eid = []
            for parent in parents:
                try:
                    parents_eid.append(node_to_C4R[parent])
                except KeyError:
                    # parent is in pedigree but not assigned G4RD ID
                    pass
            members[member]["parents_eid"] = parents_eid
        else:
            members[member]["parents_eid"] = None

    # get proband id
    proband_id = ped['pedigree']['proband']
    proband_id = node_to_C4R[proband_id]
    print(f"Proband in pedigree is {proband_id}")
    
    return members, proband_id

def get_sex(pid: int) -> str:
    """Get sex of an individual given Phenotips ID of that individual"""
    print(f"Querying Phenotips endpoint for participant {pid}")
    pid_response = requests.get(
        f"{base_url}/rest/patients/{pid}",
        auth=auth,
    )
    if pid_response.status_code == 200:
        sex = pid_response.json()["sex"]
        print(f"Query successful; participant {pid} has sex {sex}")
    else:
        # if any participant is not present in Phenotips, exit with error
        print(
            f"Query unsuccessful; participant {pid} is not present in Phenotips"
        )

    return sex

def get_affected(pid: int, members: dict) -> str:
    """Get affected status of individual"""
    affected = None
    for member in members:
        properties = member.get("properties", None)
        if properties:
            if properties.get("id", None) == pid:
                pedigree_properties = member.get("pedigreeProperties", None)
                if pedigree_properties:
                    affected = member["pedigreeProperties"].get("carrierStatus", None)
                    if not affected:
                        affected = member["pedigreeProperties"].get("carrierStatus", None)
                    break

    return affected

def get_parents(node_id: int, relationships: dict) -> list:
    """
    Given node ID, retrieve IDs of parent nodes if they exist
    """
    for relationship in relationships:
        parents = None
        if node_id in [id["id"] for id in relationship["children"]]:
            parents = [id for id in relationship["members"]]
            break

    return parents

def write_pedigree(members: dict, family: str) -> None:
    """
    Write a pedigree text file given dictionary derived from Phenotips pedigree JSON
    """
    family = family.replace("_", "")
    with open(f"/hpf/largeprojects/tgnode/sandbox/mcouse_analysis/pedigrees/DECODER/{family}_pedigree.ped", "w") as f:
        for member in members:
            family_id = family
            sample_id = member
            member = members[member]
            if member.get("sex") == "M":
                sex = 1
            elif member.get("sex") == "F":
                sex = 2
            else:
                sex = "other"
            affected = member.get("affected")
            phenotype = 2 if affected == "affected" else 0
            paternal_id = 0
            maternal_id = 0
            if member.get("parents_eid"):
                for p in member["parents_eid"]:
                    parent_sex = members[p]["sex"]
                    if parent_sex == "M":
                        paternal_id = p
                    else:
                        maternal_id = p
            # convert IDs to crg2-pacbio compatible IDs
            sample_id = sample_id.replace("_", "").replace(".", "_")
            try:
                paternal_id = paternal_id.replace("_", "").replace(".", "_")
            except:
                pass
            try:
                maternal_id = maternal_id.replace("_", "").replace(".", "_")
            except:
                pass
            f.write(
                f"{family_id} {sample_id} {paternal_id} {maternal_id} {sex} {phenotype}\n"
            )

def get_HPO_IDs(proband_id: str) -> pd.DataFrame:
        """Query G4RD phenotips to get HPO terms for the proband"""
        hpo = requests.get(
            f"{base_url}/rest/patients/{proband_id}/suggested-gene-panels",
            auth=auth,
        )

        hpo = hpo.json()
        hpo_id = []

        for row in hpo["rows"]:
            terms = row["terms"]
            for term in terms:
                hpo_id.append(term["id"])
        
        hpo_id = list(set(hpo_id))

        return hpo_id

def hpo_to_gene_mapping(hpo_ids: list) -> pd.DataFrame:
    """Map HPO terms to genes"""
    hpo_mapping = pd.read_csv("/hpf/largeprojects/tgnode/sandbox/mcouse_analysis/HPO/download/2025-05-23_genes_to_phenotype.txt", sep="\t").drop(columns=["ncbi_gene_id", "frequency", "disease_id"])
    hpo_mapping = hpo_mapping.drop_duplicates()
    hpo_agg = hpo_mapping[hpo_mapping["hpo_id"].isin(hpo_ids)].groupby("gene_symbol").agg(lambda x:  ", ".join(x)).reset_index() # get genes associated with patient HPO terms 

    return hpo_agg

def get_ensembl_from_hgnc(hpo_agg: pd.DataFrame) -> pd.DataFrame:
    """Get Ensembl IDs for genes from HGNC mapping"""
    hgnc_mapping = pd.read_csv("/hpf/largeprojects/ccmbio/ccmmarvin_shared/pacbio_longread_pilot_phase_1/annotate_SV/HPO/HGNC_ensembl_map.csv")
    hpo_agg_ens = hpo_agg.merge(hgnc_mapping, left_on="gene_symbol", right_on="hgnc_symbol", how="left")
    
    return hpo_agg_ens

def query_ensembl_batch(gene_symbols: list) -> dict:
    """Query Ensembl REST API to get gene information for multiple genes at once"""
    ensembl_url = "http://rest.ensembl.org/lookup/symbol/homo_sapiens"
    headers = {
        "Content-Type": "application/json",
        "Accept": "application/json"
    }
    payload = {"symbols": gene_symbols}
    
    try:
        response = requests.post(ensembl_url, headers=headers, json=payload)
        response.raise_for_status()
        data = response.json()
        
        # Create a mapping of gene symbol to Ensembl ID
        result = {}
        for symbol in gene_symbols:
            if symbol in data:
                result[symbol] = data[symbol]['id']
            else:
                result[symbol] = None
        return result
        
    except (requests.RequestException, JSONDecodeError) as e:
        print(f"Error querying Ensembl API: {e}")
        return {}
    
def get_dup_genes(hpo_agg_ens: pd.DataFrame) -> list:
    """Get duplicate genes from HGNC mapping"""
    genes = []
    dup_genes = []
    for gene in hpo_agg_ens["gene_symbol"].values:
        if gene not in genes:
            genes.append(gene)
        else:
            if gene not in dup_genes:
                dup_genes.append(gene)

    return dup_genes

def get_HPO_gene_mapping(hpo_ids: list) -> pd.DataFrame:
    """Get HPO gene mapping based on patient HPO terms"""
    hpo_agg = hpo_to_gene_mapping(hpo_ids) # map HPO terms to associated genes
    hpo_agg_ens = get_ensembl_from_hgnc(hpo_agg) # get Ensembl IDs for genes from HGNC mapping
    ens_ids = query_ensembl_batch(hpo_agg_ens[hpo_agg_ens["ensembl_gene_id"].isna()]["gene_symbol"].values.tolist()) # get Ensembl IDs for genes without an Ensembl ID after HGNC mapping
    hpo_agg_ens.loc[hpo_agg_ens["ensembl_gene_id"].isna(), "ensembl_gene_id"] = hpo_agg_ens.loc[hpo_agg_ens["ensembl_gene_id"].isna(), "gene_symbol"].map(ens_ids)
    dup_genes = get_dup_genes(hpo_agg_ens) # some genes are associated with multiple Ensembl IDs (possibly due to different versions over time). Query the REST API to get the latest Ensembl ID for these genes. 
    dup_ens = query_ensembl_batch(dup_genes)
    # Only update ensembl_gene_id for genes that exist in dup_ens dictionary
    mask = hpo_agg_ens["gene_symbol"].isin(dup_ens.keys())
    hpo_agg_ens.loc[mask, "ensembl_gene_id"] = hpo_agg_ens.loc[mask, "gene_symbol"].map(dup_ens)
    hpo_agg_ens = hpo_agg_ens.drop_duplicates(subset=["gene_symbol", "hpo_id", "hpo_name", "hgnc_symbol", "ensembl_gene_id"])
    # get number of patient HPO terms per gene
    hpo_agg_ens["Number of occurrences"] = hpo_agg_ens["hpo_id"].str.split(",").str.len()
    # rename columns 
    hpo_agg_ens = hpo_agg_ens.rename(columns={"hpo_id": "HPO IDs", "hpo_name": "Features",  "ensembl_gene_id": "Gene ID", "hgnc_symbol": "Gene Symbol"})
    hpo_agg_ens = hpo_agg_ens[["Gene Symbol", "Gene ID", "Number of occurrences", "Features", "HPO IDs"]]
    
    return hpo_agg_ens

def process_sample(id, fam, auth, pid_url, pedigree_url):
    response = requests.get(f"{pid_url}/{id}", auth=auth)
    try:
        pid = response.json().get('id')
        params={"action": "familyinfo", "document_id": pid}
        response = requests.get(pedigree_url, params=params, auth=auth)
        ped_json = response.json()
        members, proband_id = get_pedigree_info(ped_json)
        write_pedigree(members, fam)
        pid_proband = requests.get(f"{pid_url}/{proband_id}", auth=auth).json().get('id')
        print(pid_proband)
        HPO_ids = get_HPO_IDs(pid_proband)
        print(len(HPO_ids))
        HPO_df = get_HPO_gene_mapping(HPO_ids)
        today = date.today()
        today = today.strftime("%Y-%m-%d")
        fam = fam.replace("_", "")
        HPO_df.to_csv(f"/hpf/largeprojects/tgnode/sandbox/mcouse_analysis/HPO/DECODER/{fam}_HPO_{today}.txt", sep="\t", index=False)
    except JSONDecodeError:
        print(f"Error: did not retrieve HPO and pedigree information for {id}")

parser = argparse.ArgumentParser(description='Process sample sheet and credentials files')
parser.add_argument('-sample_sheet', help='Tab-separated sample sheet file', required=False)
parser.add_argument('-credentials', help='Credentials file containing username and password', required=True)
parser.add_argument('-sample_id', help='Single sample ID to process', required=False)
args = parser.parse_args()

credentials = pd.read_csv(args.credentials) 
username = credentials["username"][0]
password = credentials["password"][0]
auth = (username, password)
base_url = "https://genomeclinic.ccm.sickkids.ca/"
pid_url="https://genomeclinic.ccm.sickkids.ca/rest/patients/eid/"
pedigree_url=f"https://genomeclinic.ccm.sickkids.ca/get/PhenoTips/FamilyPedigreeInterface"

def main():
    if args.sample_sheet:
        sample_sheet = pd.read_csv(args.sample_sheet, sep="\t")
        sample_sheet["DECODER_family"] = sample_sheet["Decoder_ID"].str.split('.').str[0]
        for id in sample_sheet["Decoder_ID"].values:
            print(id)
            fam = id.split(".")[0]
            if '.03' in id: # proband ID
                process_sample(id, fam, auth, pid_url, pedigree_url)
    elif args.sample_id: 
        id = args.sample_id
        fam = id.split(".")[0]
        process_sample(id, fam, auth, pid_url, pedigree_url)
    else:
        print("Error: no sample ID or sample sheet provided")


if __name__ == "__main__":
    main()