import pandas as pd
import warnings
warnings.filterwarnings("ignore")
df = pd.read_csv ('Output4.csv')
table = pd.pivot_table(df, values=['gfp_mcherry','gfp_nir'], index=['well'], columns=['para'])
table.columns = [f'{j}_{i}' for i, j in table.columns]
csv = pd.read_csv('210901_flow_results.csv')
newdf = table.merge(csv, how='outer', on='well')
para_ls1 = [1,2,4,6,8,10,12,14,16,18,20]
para_ls2 = [1,2,4,6,8,10,12,14,16,18,20]

for para1 in para_ls1:
    for para2 in para_ls2:
        s = str(para1)+'_'+str(para2)+'_gfp_mcherry'
        newdf[s] = abs(newdf[s]-newdf['pct_grn_gvn_red'])

for para1 in para_ls1:
    for para2 in para_ls2:
        s = str(para1)+'_'+str(para2)+'_gfp_nir'
        newdf[s] = abs(newdf[s]-newdf['pct_grn_gvn_nir'])

min_error_mcherry = 1000000000
min_error_nir = 1000000000
min_error_mcherry_id = ''
min_error_nir_id = ''

for para1 in para_ls1:
    for para2 in para_ls2:
        s = str(para1) + '_' + str(para2) + '_gfp_mcherry'
        if newdf[s].mean() < min_error_mcherry:
            min_error_mcherry = newdf[s].mean()
            min_error_mcherry_id = s

for para1 in para_ls1:
    for para2 in para_ls2:
        s = str(para1) + '_' + str(para2) + '_gfp_nir'
        if newdf[s].mean() < min_error_nir:
            min_error_nir = newdf[s].mean()
            min_error_nir_id = s

print('min_error_mcherry')
print(min_error_mcherry)

print('min_error_nir')
print(min_error_nir)

print('min_error_mcherry_id')
print(min_error_mcherry_id)

print('min_error_nir_id')
print(min_error_nir_id)