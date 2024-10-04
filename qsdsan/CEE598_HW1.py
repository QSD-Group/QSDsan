import pandas as pd

# File path setup
faostat_file_1 = 'FAO1.xlsx'
faostat_file_2 = 'FAO2.xlsx'
faostat_file_3 = 'FAO3.xlsx'
faostat_file_4 = 'FAO4.xlsx'
faostat_file_5 = 'FAO5.xlsx'
faostat_file_6 = 'FAO6.xlsx'
faostat_file_7 = 'FAO7.xlsx'
faostat_file_8 = 'FAO8.xlsx'
faostat_file_9 = 'FAO9.xlsx'
vwc_file_path = 'VWC.xlsx'

# Loading FAOSTAT data
faostat_data_1 = pd.read_excel(faostat_file_1)
faostat_data_2 = pd.read_excel(faostat_file_2)
faostat_data_3 = pd.read_excel(faostat_file_3)
faostat_data_4 = pd.read_excel(faostat_file_4)
faostat_data_5 = pd.read_excel(faostat_file_5)
faostat_data_6 = pd.read_excel(faostat_file_6)
faostat_data_7 = pd.read_excel(faostat_file_7)
faostat_data_8 = pd.read_excel(faostat_file_8)
faostat_data_9 = pd.read_excel(faostat_file_9)
vwc_data = pd.read_excel(vwc_file_path)

# Combine all FAOSTAT data into one DataFrame
faostat_data_combined = pd.concat([faostat_data_1, faostat_data_2, faostat_data_3, 
                                   faostat_data_4, faostat_data_5, faostat_data_6,
                                   faostat_data_7, faostat_data_8, faostat_data_9], ignore_index=True)

# Cleaning up the VWC data
vwc_data_cleaned = vwc_data.rename(columns={'Product description (HS)': 'Item', 'VWC': 'VWC'})

# Merging FAOSTAT with VWC data
merged_data = pd.merge(faostat_data_combined, vwc_data_cleaned[['Item', 'VWC']], on='Item', how='inner')

# Calculating Virtual Water Trade (VWT) (Export/Import volume * VWC)
merged_data['VWT'] = merged_data['Value'] * merged_data['VWC']

# Save the merged data to a CSV file
merged_data.to_csv('output_combined_vwt_data_corrected.csv', index=False)

# Creating and saving the pivot table (using Exporting and Importing Countries)
export_import_matrix = pd.pivot_table(merged_data, values='VWT', 
                                      index='Exporting', 
                                      columns='Importing', 
                                      aggfunc='sum', fill_value=0)

# Save the pivot table to a CSV file
export_import_matrix.to_csv('export_import_matrix.csv')

# Print the first few rows of the pivot table
print(export_import_matrix.head())