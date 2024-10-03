# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 21:03:25 2024

@author: Junhyung Park
"""

import pandas as pd

# 파일 경로 설정
faostat_file_1 = 'FAOSTAT1.xlsx'
faostat_file_2 = 'FAOSTAT2.xlsx'
faostat_file_3 = 'FAOSTAT3.xlsx'
vwc_file_path = 'VWC.xlsx'

# FAOSTAT 데이터 불러오기
faostat_data_1 = pd.read_excel(faostat_file_1)
faostat_data_2 = pd.read_excel(faostat_file_2)
faostat_data_3 = pd.read_excel(faostat_file_3)

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
