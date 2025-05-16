#!/bin/bash

# Comprobación del argumento
if [ -z "$1" ]; then
	echo "Error: no ha proporcionado el directorio de entrada"
	echo "Uso: $0 <directorio_con_fasta>"
	exit 1
fi

# Variables
filedir=$1
out="${filedir}_output"

# Mensajes iniciales para el usuario
echo ""
echo -e "Directorio de entrada (alineamientos): $filedir"
echo -e "Directorio de salida (archivos generados): $out"
echo "Se generará un nuevo directorio en la ruta actual: ./$out"
echo ""
read -p "Pulse ENTER para comenzar el análisis..." intro

# Creamos el directorio de salida
mkdir -p "./$out"

# Testeo de los modelos de sustitución para cada alineamiento con iqtree2
echo -e "\nSe va a testar el mejor modelo de sustitución para cada alineamiento...\n"

for file in "./$filedir"/*; do
	base=$(basename "$file" .fasta)
    	outprefix="$out/${base}_model"

	# Si ya existe el árbol con modelo, lo saltamos
	if [[ -f "${outprefix}.treefile" ]]; then
		echo "Ya existe ${outprefix}.treefile, se omite"
		continue
	fi

	iqtree2 -s "$file" -m TEST -pre "$outprefix" -nt AUTO --keep-ident
done

# Extraemos el mejor modelo para cada alineamiento y evaluamos la señal filogenética
echo -e "\nSe va a realizar el análisis de Likelihood-mapping...\n"

for file in "./$out/"*_model.iqtree; do
	base=$(basename $file _model.iqtree)
	fasta="${filedir}/${base}.fasta"
	limap_prefix="$out/${base}_limap"

	if [[ -f "${limap_prefix}.treefile" ]]; then
		echo "Ya existe ${limap_prefix}.treefile, se omite"
        	continue
	fi

	# Extraemos el mejor modelo de sustitución
	model=$(grep "Best-fit model according to BIC:" $file | awk -F ": " '{print $2}') 
	echo "$model" > "$out/${base}.model.txt"
	read -p "Modelo óptimo para $base: $model (pulse ENTER para continuar) " enter

	# Likelihood-mapping analysis
	iqtree2 -s $fasta -lmap ALL -n 0 -m $model -pre "$limap_prefix" --keep-ident
done

# Test de topología de los árboles
echo -e "\nComparando la topología de los árboles (concatenados+gen vs. gen individual)...\n"
core_tree="$out/core_model.treefile"

for file in "./$out/"Rv*_model.treefile; do
	# Concatenación de los archivos .treefile
	base=$(basename "$file" _model.treefile)
	model=$(cat "$out/${base}.model.txt")
	cat "$file" "$core_tree" > "$out/${
base}+core.treefile"
	echo "Se ha concatenado los árboles ${base}_model.treefile y core_model.treefile"

	# Comparación por pares de la topología
	iqtree2 -s "./$filedir/${base}.fasta" -z "$out/${base}+core.treefile" -n 0 -m "$model" -zb 10000 -zw -au -pre "$out/${base}_topotest" --keep-ident -redo
done

#Eliminamos los archivos temporales
rm "$out"/*.model.txt

# Visualizamos el resultado
for file in $out/Rv*_topotest.iqtree; do
	cat $file
	base=$(basename $file _topotest.iqtree)
	read -p "Está visualizando el archivo .iqtree final de $base. Pulse ENTER para avanzar al siguiente."
done

# Redistribuimos los archivos en carpetas para su uso
lm_dir="likelihood_map"
modeltrees="model_treefile"
topo_t="topology_iqtree"

mkdir -p "$modeltrees"
mv $out/*_model.treefile "$modeltrees"

mkdir -p "$lm_dir"
mv $out/*_limap.lmap.svg "$lm_dir"

mkdir -p "$topo_t"
mv $out/*_topotest.iqtree "$topo_t"
