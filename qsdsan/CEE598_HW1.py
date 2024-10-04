# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 21:03:25 2024

@author: Junhyung Park
"""

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

# VWC 데이터 불러오기
vwc_data = pd.read_excel(vwc_file_path)

# 필요한 열만 사용하여 VWC 데이터 정리
vwc_data_cleaned = vwc_data.rename(columns={'Product description (HS)': 'Item', 'vwc': 'VWC'})

# FAOSTAT 데이터를 모두 하나로 결합
faostat_data_combined = pd.concat([faostat_data_1, faostat_data_2, faostat_data_3, 
                                   faostat_data_4, faostat_data_5, faostat_data_6,
                                   faostat_data_7, faostat_data_8, faostat_data_9], ignore_index=True)

# FAOSTAT와 VWC 데이터 매칭
merged_data = pd.merge(faostat_data_combined, vwc_data_cleaned[['Item', 'VWC']], on='Item', how='inner')

# 'Reporter Countries'와 'Partner Countries' 열 삭제
if 'Reporter Countries' in merged_data.columns:
    merged_data.drop(columns=['Reporter Countries', 'Partner Countries'], inplace=True)

# 가상 물 무역(VWT) 계산 (수출/수입량 * VWC)
merged_data['VWT'] = merged_data['Value'] * merged_data['VWC']

# 최종 데이터 출력
print(merged_data.head())

# 최종 데이터 CSV 파일로 저장하기 (원하는 경우)
merged_data.to_csv('output_combined_vwt_data.csv', index=False)

# 피벗 테이블 생성 및 저장
export_import_matrix = pd.pivot_table(merged_data, values='VWT', 
                                      index='Exporting', 
                                      columns='Importing', 
                                      aggfunc='sum', fill_value=0)

# 피벗 테이블 CSV 파일로 저장
export_import_matrix.to_csv('export_import_matrix.csv')

# 피벗 테이블의 첫 몇 행 출력
print(export_import_matrix.head())