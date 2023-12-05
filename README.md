# Diseño y validacion in silico de cebadores especificos para peces

_Martin Holguin Osorio_\
_Enero de 2020_ 

>Estos experimentos estan basados en el trabajo hecho por Eric Coissac en 2015:\
>[Designing a new DNA metabarcode for fish](https://metabarcoding.org/IMG/html/primerdesign.html)

## Dependencias

Para correr y analizar los datos descargados se requieren los siguientes programas:\
* Un ordenador con un sistema operativo UNIX, o un servidor o maquina virtual que permita ejecutar codigo en linea de comandos.
* El ordenador o maquina debe tener instalados los siguientes programas:
    * Paquete de programas [OBITools](https://pythonhosted.org/OBITools/welcome.html#installing-the-obitools).
    * Programa [ecoPCR](https://git.metabarcoding.org/obitools/ecopcr/-/wikis/home)
    * Programa [ecocebadores](https://git.metabarcoding.org/obitools/ecocebadores/-/wikis/home)
    * Programa R [r-project](https://www.r-project.org/) y tener los paquetes ROBITools, ROBITaxonomy y ROBIBarcodes.

Una vez este todo listo, iniciamos:

###Descarga y formateo de los datos

```
# Creo carpeta para almacenar todos los datos desde mi carpeta personal
mkdir mitocondria

# entro a la carpeta
cd mitochondria

# descargo secuencias de DNA mitcondrial completas del ncbi
wget "ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/mitochondrion/mitochondrion.1.genomic.gbff.gz"

wget "ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/mitochondrion/mitochondrion.2.genomic.gbff.gz"

#descargo taxonomia asociada a estas secuencias
wget 'ftp://ftp.ncbi.nlm.nih.gov://pub/taxonomy/taxdump.tar.gz'

#extraccion de la taxonomia en una carpeta 
mkdir ncbi20150518
cd ncbi20150518/
tar xf ../taxdump.tar.gz
cd ..
```

### Formateo de los datos para OBITools

```
#activo a obitools desde donde lo tenga instalado
cd /home/usuario/programas
./obitools

#vuelvo a la direccion de trabajo
cd /home/usuario/mitochondria2

#formateo los datos para obitools
obitaxonomy -t ncbi20150518 -d ncbi20150518

#uno todos los genomas en un archivo fasta unico y verifico
obiconvert mitochondria/* > mito.all.fasta
head -5 mito.all.fasta

#busco los vertebrados en los archivos de taxonomia en base a su taxid
ecofind -d ncbi20150518 '^vertebrata$'

#miro errores existentes "como el de un genero llamado vertebrata"
ecofind -d ncbi20150518 -p 1261581

#corrigo esos errores "reanoto y selecciono los genomas"
obiannotate -d ncbi20150518 \
            --with-taxon-at-rank=species \
            mito.all.fasta | \
obiannotate -S 'ori_taxid=taxid' | \
obiannotate -S 'taxid=species' | \
obiuniq -c taxid | \
obiselect -c taxid -n 1 -f count -M > mito.one.fasta

#verifico numero de secuencias tras filtrado
obicount mito.all.fasta

#selecciono los genomas de vertebrata
obigrep -d ncbi20150518 -r 7742 mito.one.fasta > mito.vert.fasta

#formateo archivo fasta en la base de datos de ecoPCR
obiconvert -d ncbi20150518 --ecopcrdb-output=mito.vert mito.vert.fasta

#miro a teleostei en base a su taxid
ecofind -d ncbi20150518 '^Teleostei$'
```

### Estimacion de las mejores parejas de cebadores para el grupo de interes (Teleostei)

```
#estimo los mejores sets de cebadores con ecocebadores y selecciono el mejor par que me da el archivo generado

#desde la carpeta src del programa "en mi caso es /home/usuario/programas/ecopcr"
./ecocebadores -d /home/usuario/mitochondria2/mito.vert  -e 3 -3 2  -l 30 -L 150 -r 32443 -c > Teleostei.ecocebadores
```

### Probamos la nueva pareja de cebadores

```
#simulo la PCR con los dos mejores cebadores dados en el archivo Teleostei.ecocebadores 

#desde la carpeta src del programa "en mi caso es /home/usuario/programas/ecocebadores/src"
./ecoPCR -d /home/usuario/mitochondria2/mito.vert -e 5 -l 30 -L 300 -c ACACCGCCCGTCACTCTC ACCTTCCGGTACACTTAC > Teleostei.ecoPCR
```

### Iniciamos sesion en R


```

#Abro R y cargo paquetes
R
library(ROBITools)
library(ROBITaxonomy)
library(ROBIBarcodes)

#cargo los datos resultantes de ecoPCR junto con su taxonomia
fish = read.ecopcr.result('Teleostei.ecoPCR')
taxo = read.taxonomy('ncbi20150518')
#verifico
head(fish,n = 2)

#identifico que secuencias son de peces en base a los datos del archivo de taxonomia
teleo.taxid = ecofind(taxo,'^Teleostei$')
teleo.taxid

#identifico que secuencias son de peces entre los resultados de ecoPCR
is_a_fish=is.subcladeof(taxo,fish$taxid,teleo.taxid)
table(is_a_fish)

#Pruebo la conservación de los sitios de ligamiento de los cebadores
Fish.forward = ecopcr.forward.shanon(ecopcr = fish,
                                     group = is_a_fish)
Fish.reverse = ecopcr.reverse.shanon(ecopcr = fish,
                                     group = is_a_fish)
									 
#grafico los resultados
pdf("cebadores_Teleoistei.pdf")
par(mfcol=c(2,2))
dnalogoplot(Fish.forward$'TRUE',
            primer = "ACACCGCCCGTCACTCTC",
            main='Forward Fish')
dnalogoplot(Fish.forward$'FALSE',
            primer = "ACACCGCCCGTCACTCTC",
            main='Forward not Fish')

dnalogoplot(Fish.reverse$'TRUE',
            primer = "ACCTTCCGGTACACTTAC",
            main='Reverse Fish')
dnalogoplot(Fish.reverse$'FALSE',
            primer = "ACCTTCCGGTACACTTAC",
            main='Reverse not Fish')
dev.off()			

			
#grafico desaciertos en los cebadores
pdf("mismatches.cebadores_Teleostei.pdf")
par(mfcol=c(1,1))
mismatchplot(fish,group = is_a_fish,
             legend=c('2722 vertebrados diferentes a peces','2617 peces'))+title(xlab="Numero de desaciertos en el primer forward", 
                                   ylab="Numero de desaciertos en el primer reverse",
                                   main = 'Distribucion del numero de desaciertos en la pareja de cebadores')
dev.off()			 
			 
#miro cual es la resolucion taxonomica de los cebadores	
only.fish=fish[is_a_fish,]

res = resolution(taxo,only.fish)
resolution = with(only.fish,
                  unique(data.frame(species_name,taxid,rank=res))
                 )
t(t(sort(table(resolution$rank)/length(resolution$rank),decreasing = TRUE)))


##Mismo analisis pero para familias

#verifico amplificaciones en Cyprinidae
#creo un grupo que contenga solo a Cyprinidae (en vez de teleostei) en base a su id
is_a_cyprinidae = is.subcladeof(taxo,fish$taxid,7953)
only.cyprinidae=fish[is_a_cyprinidae,]

#Verifico resolucion de los cebadores entre Cyprinidae (en vez de teleostei)
res = resolution(taxo,only.cyprinidae)
resolution = with(only.cyprinidae,
                  unique(data.frame(species_name,taxid,rank=res))
                 )
t(t(sort(table(resolution$rank)/length(resolution$rank),decreasing = TRUE)))

#verifico resolucon de cebadores para los otros peces afuera de Cyprinidae
no.cyprinidae=fish[is_a_fish & !is_a_cyprinidae,]
res = resolution(taxo,no.cyprinidae)
resolution = with(no.cyprinidae,
                  unique(data.frame(species_name,taxid,rank=res))
                 )
t(t(sort(table(resolution$rank)/length(resolution$rank),decreasing = TRUE)))
```
