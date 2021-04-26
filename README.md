# gTidesKozai
#Compilation:
gcc gTidesKozai.c -m -lgsl -lgslcblas -o gTides

#Running the code:
./gTides FILENAME

#Plotting Results
python3 plotcmp.py FILENAME
