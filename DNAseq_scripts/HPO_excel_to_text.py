import pandas as pd

# Can't export HPO from Epic directly for clinical samples, so HPO terms have been saved as sheets in an excel file

excel_file = pd.ExcelFile('/hpf/largeprojects/tgnode/sandbox/mcouse_analysis/HPO/genesteps/genesteps_HPO_01_to_50.xlsx')

# Get sheet names
sheet_names = excel_file.sheet_names

# Create a dictionary to store each sheet's data
sheets_data = {}

# Loop through each sheet and store its data in the dictionary
for sheet_name in excel_file.sheet_names:
    # Read the sheet into a DataFrame
    df = excel_file.parse(sheet_name, names=["Gene Symbol", "Gene ID", "Number of occurrences", "Features", "HPO IDs"] , header=None)
    if len(df) != 0:
        # Filter out rows where any column starts with '#'
        df = df[~df.apply(lambda row: row.astype(str).str.startswith('#').any(), axis=1)]
        df = df[~df["Gene Symbol"].isna()]

        # Store the cleaned DataFrame in the dictionary
        sheets_data[sheet_name] = df
# Now you have each sheet's data stored in the sheets_data dictionary
# Access each sheet's data using sheets_data['sheet_name']
for sheet in sheets_data:
    fam = "RGSE" + sheet
    sheets_data[sheet].to_csv(f"/hpf/largeprojects/tgnode/sandbox/mcouse_analysis/HPO/genesteps/{fam}_HPO.txt", sep="\t", header=None, index=False)
