#!/bin/zsh
# Download the public proteomics cohorts that need no login into data/proteomics/<study>/raw and
# write data/proteomics/downloads.tsv (study, file, url, bytes, checksum type, expected, observed, match, sha256, time of the check; files present are verified, not downloaded again).
# Vijayakumar2022 (Dryad) is not here: the Dryad file API needs a token and the landing page refuses scripted access.
set -u
cd "$(dirname $0)/../.."
LOG=data/proteomics/downloads.tsv
mkdir -p data/proteomics
[[ -f $LOG ]] || printf "study\tfile\turl\tbytes\tchecksum_type\texpected\tobserved\tmatch\tsha256\tverified\n" > $LOG
get() {  # study file url checksum_type expected; a checksum mismatch is retried (3 attempts)
  local d=data/proteomics/$1/raw; mkdir -p "$d"
  local f="$d/$2" obs sha m try
  for try in 1 2 3; do
    if [[ ! -s "$f" ]]; then curl -sfL --retry 3 -o "$f.part" "$3" && mv "$f.part" "$f" || { rm -f "$f.part"; echo "FAILED $1 $2 (attempt $try)"; sleep 20; continue; }; fi
    sha=$(shasum -a 256 "$f" | cut -d' ' -f1)
    case $4 in
      md5) obs=$(md5 -q "$f");;
      sha1) obs=$(shasum -a 1 "$f" | cut -d' ' -f1);;
      sha256) obs=$sha;;
      *) obs=;;
    esac
    m=n.a.; [[ -n $5 ]] && { [[ $obs == $5 ]] && m=yes || m=NO; }
    [[ $m != NO ]] && break
    echo "checksum mismatch $1 $2 (attempt $try)"; rm -f "$f"; sleep 20
  done
  [[ -s "$f" ]] || return 1
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$1" "$2" "$3" "$(stat -f %z "$f")" "${4:-none}" "${5:-}" "$obs" "$m" "$sha" "$(date '+%Y-%m-%d %H:%M:%S')" >> $LOG
}
Z=https://zenodo.org/api/records/14772494/files
get Buggert2025 meta_UKsamples.csv $Z/meta_UKsamples.csv/content md5 a348aeeaaae4e2ad7d068fb3b6472a04
get Buggert2025 meta_SEsamples.csv $Z/meta_SEsamples.csv/content md5 0e005ba4c16781c8da22d26a11b1e7d2
get Buggert2025 plate1_olink.csv $Z/plate1_olink.csv/content md5 8c07ec6cf0a40497785d6c887d1fdac8
get Buggert2025 plate2_olink.csv $Z/plate2_olink.csv/content md5 6d3f00ef9cf1a51fa0f3d1c273e6ffa9
get Buggert2025 plate3_olink.csv $Z/plate3_olink.csv/content md5 e3cbb63d6e875f087d0661f30902f24c
M=https://data.mendeley.com/public-files/datasets
get Cervia2024 Cervia_Hasler_et_al_2024_Proteomics_Clinical_Data_6M_timepoint.xlsx $M/dvf6yvrg4x/files/9a3eedfb-93b8-4b43-8dbb-ffe8adc45b1e/file_downloaded sha256 3752cf6bf441e0781de5717f76d8f6abd3d5ecc2a5dd89d9c8624ae2aafd85ea
get Cervia2024 Cervia_Hasler_et_al_2024_Proteomics_Clinical_Data_Acute_plus_healthy-81EAkK.xlsx $M/dvf6yvrg4x/files/ec87a7da-b850-4028-b024-053f78445398/file_downloaded sha256 c8f3226b9c1c2be744ed4b65843a6e6f1c876d405f1966888ba5eb3e20c7364b
get Cervia2024 Cervia_Hasler_et_al_2024_Proteomics_Clinical_Data_Labels-TAklSO.xlsx $M/dvf6yvrg4x/files/859c49d0-ae4c-4fc1-85c2-8d178e274fa3/file_downloaded sha256 6edeea54657f524d2dd9cbb9ede9c49fc9a1d030f4cdaafcd29dd914dde58772
get Cervia2024 Cervia_Hasler_et_al_2024_Proteomics_SomaScan_Aptamer_Info-KnqtVJ.xlsx $M/dvf6yvrg4x/files/6284c5f0-33b7-457a-8bf2-642f46d69184/file_downloaded sha256 68ae0202f2b75c5c37e37f42ee888edf39ca305bf1b2a81e600dc923e8120f2b
get Cervia2024 Cervia_Hasler_et_al_2024_Mass_Spectrometry_6M_timepoint-VOZ1Hu.xlsx $M/dvf6yvrg4x/files/530d6021-7fe1-4c95-8d09-f525bba17c73/file_downloaded sha256 cf62ee6ca935b4e8804b1a6b7a7aac408cdbf1d59d4b8bf15a0dbabd1f5dc1c9
get Su2022 "Table S1.xlsx" $M/96v329bg7g/files/dbfc7e42-66cb-4529-b863-c792bbed122d/file_downloaded sha256 2a840668db7c8c6bb37aa486f24d6db2433329a068efeb8a2acd92b3a64e1a88
get Su2022 "Table S2.xlsx" $M/96v329bg7g/files/482dc674-8fa1-494d-a156-2341c43be549/file_downloaded sha256 611a40a9dd2753c81975a15ffd4e6c8f8d535a492bcefca4ac2548122ca4d828
P=https://ftp.pride.ebi.ac.uk/pride/data/archive
get Kovarik2023_plasma MQ_Search.zip $P/2023/03/PXD036969/MQ_Search.zip sha1 b4b7b8950304ac704bea52bb0acc8f4e778f40b4
get Kovarik2023_macrophage MQ_Search.zip $P/2023/03/PXD036970/MQ_Search.zip sha1 f2efe36f8ad924270cc2e0e119c1f22bad5dc4f2
get Captur2022 Supplementary_Data_1.xlsx $P/2022/10/PXD036590/Supplementary_Data_1.xlsx sha1 ffa301f710c6163fe42b2c121d6c55568554cfec
get Captur2022 Plate_1_Injection1_TransitionResults_QuantQuali_YE_IS.csv $P/2022/10/PXD036590/Plate_1_Injection1_TransitionResults_QuantQuali_YE_IS.csv sha1 55de96582d12cd2bec418f4f1ae29a738f4f114f
get Captur2022 Plate_1_Injection_2_TransitionResults_QuantQuali_YE_IS.csv $P/2022/10/PXD036590/Plate_1_Injection_2_TransitionResults_QuantQuali_YE_IS.csv sha1 ee864d2802233a527dbc4ebcc217fee54c11f7bb
get Captur2022 Plate_2_Injection1_TransitionResults_QuantQuali_YE_IS.csv $P/2022/10/PXD036590/Plate_2_Injection1_TransitionResults_QuantQuali_YE_IS.csv sha1 7d1d9fa6a7b2e9e96c8123811a8f22b42496ec87
get Captur2022 Plate_2_Injection2_TransitionResults_QuantQuali_YE_IS.csv $P/2022/10/PXD036590/Plate_2_Injection2_TransitionResults_QuantQuali_YE_IS.csv sha1 d3f8abbcd354a1036d8734753da618381a567e74
get Captur2022 Plate_3_Injection1_TransitionResults_QuantQuali_YE_IS.csv $P/2022/10/PXD036590/Plate_3_Injection1_TransitionResults_QuantQuali_YE_IS.csv sha1 7811dc82f99927425e86d25b4de61955e17834ee
get Captur2022 Plate_3_Injection2_TransitionResults_QuantQuali_YE_IS.csv $P/2022/10/PXD036590/Plate_3_Injection2_TransitionResults_QuantQuali_YE_IS.csv sha1 9c0e8adf62a6fd91a6147ecdf2f631f8bf8b543c
get Captur2022 Plate_4_Injection1_TransitionResults_QuantQuali_YE_IS.csv $P/2022/10/PXD036590/Plate_4_Injection1_TransitionResults_QuantQuali_YE_IS.csv sha1 a368f895aa31b22c366af04349932b4df39ed3b1
get Captur2022 Plate_4_Injection2_TransitionResults_QuantQuali_YE_IS.csv $P/2022/10/PXD036590/Plate_4_Injection2_TransitionResults_QuantQuali_YE_IS.csv sha1 db6f13874986a55e7fa45ee184e864b1ffc45b7a
get Chowdhury2024 metadata.sdrf.tsv $P/2024/10/PXD047623/metadata.sdrf.tsv none ""
get Chowdhury2024 PXD047623_community_annotated.sdrf.tsv $P/2024/10/PXD047623/PXD047623_community_annotated.sdrf.tsv none ""
get Chowdhury2024 diann.7z $P/2024/10/PXD047623/diann.7z none ""
# only the DIA-NN summary tables are taken out of diann.7z (it also holds a 3.5 GB spectral library);
# needs a python with py7zr (PY7Z=/path/to/python; 7-Zip is not installed here)
if [[ ! -s data/proteomics/Chowdhury2024/raw/diann/report.pg_matrix.tsv ]]; then
  ${PY7Z:-python3} -c "import py7zr; py7zr.SevenZipFile('data/proteomics/Chowdhury2024/raw/diann.7z').extract(path='data/proteomics/Chowdhury2024/raw/diann', targets=['report.pg_matrix.tsv','report.pr_matrix.tsv','report.gg_matrix.tsv','report.unique_genes_matrix.tsv','report.stats.tsv','report.log.txt'])" || echo "diann.7z not extracted: install py7zr"
fi
# Zoodsma2022: the GitHub repository at the commit recorded in docs/proteomics_inventory.csv
R=data/proteomics/Zoodsma2022/raw/Postcovid-targeted-proteomics
if [[ ! -d $R ]]; then
  git clone -q https://github.com/CiiM-Bioinformatics-group/Postcovid-targeted-proteomics.git $R && git -C $R checkout -q d45e36077d735c242788e72b871ce4a3d1e728a9
fi
printf "Zoodsma2022\t(git clone)\thttps://github.com/CiiM-Bioinformatics-group/Postcovid-targeted-proteomics\t%s\tgit_commit\td45e36077d735c242788e72b871ce4a3d1e728a9\t%s\t%s\t\t%s\n" "$(du -sk $R | cut -f1)kB" "$(git -C $R rev-parse HEAD)" "$([[ $(git -C $R rev-parse HEAD) == d45e36077d735c242788e72b871ce4a3d1e728a9 ]] && echo yes || echo NO)" "$(date '+%Y-%m-%d %H:%M:%S')" >> $LOG
echo done
