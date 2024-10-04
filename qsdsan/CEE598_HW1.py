# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 21:03:25 2024

@author: Junhyung Park
"""

import pandas as pd
import matplotlib.pyplot as plt
# 파일 경로 설정
faostat_file_1 = 'FAO1.xlsx'
faostat_file_2 = 'FAO2.xlsx'
faostat_file_3 = 'FAO3.xlsx'
faostat_file_4 = 'FAO4.xlsx'
faostat_file_5 = 'FAO5.xlsx'
faostat_file_6 = 'FAO6.xlsx'

vwc_file_path = 'VWC.xlsx'

# FAOSTAT 데이터 불러오기
faostat_data_1 = pd.read_excel(faostat_file_1)
faostat_data_2 = pd.read_excel(faostat_file_2)
faostat_data_3 = pd.read_excel(faostat_file_3)
faostat_data_4 = pd.read_excel(faostat_file_4)
faostat_data_5 = pd.read_excel(faostat_file_5)
faostat_data_6 = pd.read_excel(faostat_file_6)

# VWC 데이터 불러오기
vwc_data = pd.read_excel(vwc_file_path)

# 필요한 열만 사용하여 VWC 데이터 정리
vwc_data_cleaned = vwc_data.rename(columns={'Product description (HS)': 'Item', 'vwc': 'VWC'})

# 세 개의 FAOSTAT 데이터를 하나로 결합
faostat_data_combined = pd.concat([faostat_data_1, faostat_data_2, faostat_data_3], ignore_index=True)

# FAOSTAT와 VWC 데이터 매칭
merged_data = pd.merge(faostat_data_combined, vwc_data_cleaned[['Item', 'VWC']], on='Item', how='inner')

# 가상 물 무역(VWT) 계산 (수출/수입량 * VWC)
merged_data['VWT'] = merged_data['Value'] * merged_data['VWC']

# 최종 데이터 출력
print(merged_data.head())

# 최종 데이터 CSV 파일로 저장하기 (원하는 경우)
merged_data.to_csv('output_combined_vwt_data.csv', index=False)
#%%
# 피벗 테이블 생성 및 저장
export_import_matrix = pd.pivot_table(merged_data, values='VWT', 
                                      index='Exporting Countries', 
                                      columns='Importing Countries', 
                                      aggfunc='sum', fill_value=0)

# 피벗 테이블 CSV 파일로 저장
export_import_matrix.to_csv('export_import_matrix.csv')

# 피벗 테이블의 첫 몇 행 출력
print(export_import_matrix.head())