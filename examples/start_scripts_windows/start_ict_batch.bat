# echo Isotopologe
..\..\ict.pl -c ..\data_dir\batch\chem_batch1.csv -i ..\data_dir\batch\nat_isotope.txt -m ..\data_dir\batch\measured_Isotopologe.csv -o ..\data_dir\batch\corrected_Isotopologe.csv

# echo Isotopomere
..\..\ict.pl -c ..\data_dir\batch\chem_batch2.csv -i ..\data_dir\batch\nat_isotope.txt -m ..\data_dir\batch\measured_Isotopomere.csv -o ..\data_dir\batch\corrected_Isotopomere.csv

..\..\ict.pl -c ..\data_dir\batch\chem_batch2.csv -i ..\data_dir\batch\nat_isotope.txt -m ..\data_dir\batch\measured_CI_Flux.csv -o ..\data_dir\batch\corrected_CI_Flux.csv
