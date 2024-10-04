import pandas as pd

# 파일 경로 설정
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

# FAOSTAT 데이터 불러오기
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

# FAOSTAT 데이터를 모두 하나로 결합
faostat_data_combined = pd.concat([faostat_data_1, faostat_data_2, faostat_data_3, 
                                   faostat_data_4, faostat_data_5, faostat_data_6,
                                   faostat_data_7, faostat_data_8, faostat_data_9], ignore_index=True)

# VWC 데이터 정리
vwc_data_cleaned = vwc_data.rename(columns={'Product description (HS)': 'Item', 'VWC': 'VWC'})

# FAOSTAT와 VWC 데이터 매칭
merged_data = pd.merge(faostat_data_combined, vwc_data_cleaned[['Item', 'VWC']], on='Item', how='inner')

# 'Exporting Countries'와 'Importing Countries'가 있는지 확인하고, 없을 경우 열을 추가하거나 수정
if 'Exporting Countries' not in merged_data.columns:
    print("Warning: 'Exporting Countries' column is missing!")
if 'Importing Countries' not in merged_data.columns:
    print("Warning: 'Importing Countries' column is missing!")

# 가상 물 무역(VWT) 계산 (수출/수입량 * VWC)
merged_data['VWT'] = merged_data['Value'] * merged_data['VWC']

# 데이터 CSV 파일로 저장
merged_data.to_csv('output_combined_vwt_data_corrected.csv', index=False)

# 데이터 확인을 위해 첫 몇 행 출력
print(merged_data.head())