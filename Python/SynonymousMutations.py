import matplotlib.pyplot as plt
import openpyxl as xl
import os
%matplotlib inline

threshold = 250
threshold = str(threshold)

def OpenExcel (file_name):
    os.chdir('F:/Promotion/Matlab/Data/NGS_7.2Mix_sorted_BT_CT/2BT/Synonymous Mutations')
    wb = xl.load_workbook(file_name,data_only=True, read_only=True)
    ws = wb['Results']
    data = []
    count = 1
    for row in ws.rows:
        if count < 5:
            count += 1
            continue
        for cell in row:
            try: 
                float(cell.value)
            except:
                continue
            if not cell.value == 0:
                data.append(cell.value)

    return data

d1 = OpenExcel('2BT1_syn-mut_enrichment_'+threshold+'.xlsx')
d2 = OpenExcel('2BT2_syn-mut_enrichment_'+threshold+'.xlsx')
d3 = OpenExcel('2BT3_syn-mut_enrichment_'+threshold+'.xlsx')
d4 = OpenExcel('2BT4_syn-mut_enrichment_'+threshold+'.xlsx')
            


# notched plot
#plt.figure()
#plt.boxplot(data, notch=True, showfliers=False)



data = [d1, d2, d3, d4]
# multiple box plots on one figure
plt.figure()
plt.boxplot(data, notch=True, showfliers=False, labels=['HI','WP','SL','LO'] )
plt.xticks(size = 20)
plt.yticks(size = 16)
plt.ticklabel_format(style='sci', axis='y', scilimits=(-3,-3), useMathText=True, useLocale=True)
plt.ylabel(r'Normalised enrichment x $10^{-3}$', size = 16)

#plt.show()
plt.savefig(threshold+'.png', bbox_inches='tight')
