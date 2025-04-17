import os
import pandas as pd
import scanpy as sc
from multiprocessing import Pool
import numpy as np
import gc


wdir = 'D:/pythonProject/spinal/250416/'
os.chdir(wdir)
pwd = 'F:/bz_new_h5ad_counts/remove/'


cellmarker = {
    'Oligodendrocyte Progenitor Cells': ['OLIG1', 'OLIG2'],
    'Radial Glial Cell': ['HES1'],
    'Glial Progenitor Cells': ['SPARCL1', 'HOPX'],
    'Neural Progenitor Cells': ['PAX6'],
    'Astrocytes': ['AQP4', 'GFAP', 'SLCO1C1'],
    'Microglia': ['PTPRC', 'P2RY12', 'C1QB', 'C1QC', 'CSF1R', 'CD68', 'AIF1'],
    'Inhibitory Neuron': ["GAD1", "GAD2", "SLC6A1", "GABBR2", "SLC32A1", "DLX2", "DLX1"],
    'Excitatory Neuron': ["NEUROD2", "NEUROD6", "SLC17A7","SLC17A8", "EOMES"],
    'Oligodendrocytes': ['MBP']
}

Age_map = {'WCL6W':'CS12-13','CW5W':'CS12-13','ZMD6W':'CS14-15','SLQ7W':'CS14-15','GY7W':'CS14-15','MXY8W':'CS17','JY9W':'CS18','JXL7W':'CS16','HP8W':'CS19','LHQ8W':'CS19','LY9W':'CS20','GMJ9W':'CS20','HQM10W':'CS23','Embryo10W':'CS23'}

Brain_map = {
    'Forebrain': ['Fb', 'Or', 'Die', 'Hy', 'Pall', 'SPall', 'Pall VZ', 'SPall VZ'],
    'MidBrain': ['Mb VZ', 'Mb'],
    'HindBrain': ['Hb', 'Hb VZ', 'Cere'],
    'SpinalCord': ['SpC']
}


def process_file(file):
    if not file.endswith(".h5ad"):
        return []

    results = []
    try:
        ad = sc.read(os.path.join(pwd, file))

        sample = f"{file.split('_')[0]}_{file.split('_')[1]}_{file.split('_')[2]}"
        age = Age_map[f"{file.split('_')[0]}{file.split('_')[1]}"]
        print(sample)
        print(age)

        available_genes = set(ad.var_names)

        for region, sub_regions in Brain_map.items():
            region_mask = ad.obs['Brain Zone'].isin(sub_regions)
            region_ad = ad[region_mask]
            region_total = region_ad.n_obs

            for sub_region in sub_regions:
                if sub_region not in ad.obs['Brain Zone'].values:
                    continue

                sub_mask = ad.obs['Brain Zone'] == sub_region
                sub_ad = ad[sub_mask]
                sub_total = sub_ad.n_obs

                for cell_type, markers in cellmarker.items():
                    cell_markers = [m for m in markers if m in available_genes]
                    cell_count = 0

                    if cell_markers:
                        cell_ad = sub_ad[:, cell_markers].copy()
                        sc.pp.filter_cells(cell_ad, min_genes=1)
                        cell_count = cell_ad.n_obs
                        del cell_ad

                    sub_prop = cell_count / sub_total if sub_total else 0
                    region_prop = cell_count / region_total if region_total else 0

                    results.append([
                        sample, age, region, sub_region, f"{age}_{sub_region}",
                        cell_type, region_total, sub_total, cell_count,
                        region_prop, sub_prop
                    ])

                del sub_ad, sub_mask
                gc.collect()

            del region_ad, region_mask
            gc.collect()

        del ad
        gc.collect()

    except Exception as e:
        print(f"Error processing {file}: {str(e)}")
    finally:
        return results




def main():
    files = sorted([f for f in os.listdir(pwd) if f.endswith(".h5ad")])

    with Pool(processes=min(3, os.cpu_count())) as pool:
        results = pool.map(process_file, files)

    df = pd.DataFrame(
        [item for sublist in results for item in sublist],
        columns=['sample', 'age', 'BrainZone', 'Brainsubstructures', 'age_BrainZone',
                 'celltype', 'BrainZonesum', 'Brainsubstructuressum',
                 'cell_sum', 'cell_allpro', 'cell_pro']
    )
    del results
    gc.collect()
    lists = [['age', 'BrainZone', 'BrainSubZone', 'celltype', 'BrainZonesum', 'BrainSubZonesum', 'cell_sum', 'cell_pro',
              'sub_cell_sum', 'sub_cell_pro']]

    for celltype in list(set(df['celltype'])):
        df1 = df[df['celltype'].isin([celltype])]
        for age in list(set(df1['age'])):
            df2 = df1[df1['age'].isin([age])]
            for BrainZone in list(set(df2['BrainZone'])):
                df3 = df2[df2['BrainZone'].isin([BrainZone])]
                BrainZonesum = df3['Brainsubstructuressum'].sum()
                cell_sum = df3['cell_sum'].sum()
                cell_pro = cell_sum / BrainZonesum
                for BrainSubZone in list(set(df3['Brainsubstructures'])):
                    df4 = df3[df3['Brainsubstructures'].isin([BrainSubZone])]
                    BrainSubZonesum = df4['Brainsubstructuressum'].sum()
                    sub_cell_sum = df4['cell_sum'].sum()
                    sub_cell_pro = sub_cell_sum / BrainZonesum
                    lists.append(
                        [age, BrainZone, BrainSubZone, celltype, BrainZonesum, BrainSubZonesum, cell_sum, cell_pro,
                         sub_cell_sum, sub_cell_pro])
    a = np.asarray(lists)
    pd.DataFrame(a).to_csv('summary_cellbin_bar.csv')

    del df
    gc.collect()


if __name__ == '__main__':
    main()

