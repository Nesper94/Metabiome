.. Metabiome documentation master file, created by
   sphinx-quickstart on Tue Jan 26 23:09:41 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Metabiome's documentation!
=====================================

Metabiome es un pipeline para el estudio de comunidades microbianas derivadas de un enfoque
metagenómico. Puede ser implementado con cualquier tipo de muestra ambiental.


La metagenómica o “shotgun metagenomics” (de ahora en adelante solo metagenómica) es una técnica de
análisis de ADN no dirigida, de todos los ("meta") genomas microbianos ("genómica") presentes en una
muestra ambiental (Jansson y col., 2018).

En términos generales un estudio de metagenómica comprende los siguientes pasos: (i) muestreo,
extracción y secuenciamiento; (ii) control de calidad o pre-procesamiento de las lecturas; (iii)
ensamblaje; (iv) análisis de las secuencias ensambladas para perfilar las características
taxonómicas y funcionales del microbioma; y finalmente (v) el análisis estadístico de los datos
(Quince y col., 2017).


**1. Implementación:** Metabiome fue construido bajo la consideración de diferentes prácticas bioinformáticas para el estudio de metagenomas, contenidos en diferentes revisiones especializadas (Quince y col., 2017; Knight y col., 2018; Reiter, Brooks y Reid, 2020). Está escrito en el lenguaje de programación Bash, y puede ser instalado en cualquier distribución de GNU/Linux, y en otros sistemas operativos como Windows y macOS. Utiliza ambientes de Conda, el sistema de gestión de paquetes de código abierto y de entornos, facilitando la instalación de softwares para usuarios no-administradores. Las herramientas del procesamiento están divididas por módulos localizados en ambientes diferentes nombrados de acuerdo al tipo de proceso que se desee realizar **(insertar driagrama de flujo de metabiome)**.  Metabiome puede ser descargado desde el repositorio de github disponible en línea **(crear enlace para github)**. Los requerimientos para utilizar este pipeline son un conocimiento básico de lenguaje de programación o ejecución de órdenes mediante línea de comandos, el software Anaconda, y acceso a un clúster remoto con capacidad para ejecutar procesos muy demandantes tanto de tiempo como de memoria.

**2. Modularidad**
Se diseñaron XNÚMERO de módulos que contienen las herramientas necesarias para ejecutar los diferentes puntos en un análisis de metagenómica. Están separados mediante ambientes de conda, creados a partir de un archivo .yaml, en el que se describen los softwares que implementa cada uno y la versión requerida. Dichos archivos están almacenados en la carpeta conda_envs. Un módulo puede tener una o más de una herramienta o software, y cada una tiene un script por separado para su ejecución, y están almacenados en el directorio “scripts” **(insertar pantallazo de terminal con los ambientes .yaml)**
Para la instalación es necesario ejecutar en la terminal el archivo install.sh, que corresponde al script de instalación del pipeline, es decir, con el que se crean los ambientes de conda, y se instala de cada software dentro de su respectivo ambiente. Cada software fue seleccionado después de una larga revisión bibliográfica. Algunos de ellos fueron tenidos en cuenta gracias las consideraciones de la Evaluación Crítica de la Interpretación del Metagenoma (CAMI por sus siglas en inglés) (Sczyrba y col., 2017). Esta es una iniciativa impulsada por la comunidad que evalúa los métodos computacionales para el análisis del metagenoma de manera integral y más objetiva, en la que se incluyen conjuntos de datos de diferentes muestras ambientales y de diferente complejidad.

    **2.1 Metabiome preprocessing**
    **2.3 Metabiome genome assembly**
    **2.4 Metabiome taxonomic binnig**
    **2.5 Metabiome taxonomic profiling**
    **2.6 Metabiome humman2**
    **2.7 Metabiome picking 16S**
    



**3. Uso de Metabiome:** Para la ejecución de Metabiome es necesario escribir en una terminal el script de “corrida” denominado metabiome.sh, seguido del nombre del script que ejecuta el programa que se desee, la dirección de la carpeta input (opción -i seguida de la dirección) con las secuencias a analizar y la dirección de la carpeta output (opción -o seguida del nombre de la carpeta output) que es donde se espera que se almacenen los datos. Cada script contiene una función de ayuda que permite ver los parámetros necesarios para la ejecución del software, y algunas opciones adicionales del mismo con su respectiva descripción, brindando mayor facilidad de uso. Esta ayuda puede ser visualizadas indicando nombre del script (.sh) que corre dicho software, seguido del comando -h (ó –h) **(insertar pantallazo con la función de ayuda y con ejemplo de cómo correr un script)**, y las opciones empleadas mediante la opción -opt, que son parámetro opcional. Así, el usuario puede ejecutar cada software desde la visualización inicial de las secuencias y seguir el flujo de trabajo y, por otro lado, si por ejemplo tiene secuencias ya depuradas, puede utilizarlas para correr cualquiera de las herramientas de análisis posterior al preprocesamiento. Además, da la posibilidad de emplear otras características y opciones de cada herramienta sin modificar el código. 


**4. Archivos de entrada o input:**
Los archivos de entrada consiste en un paquete de datos que puede tener el formato .fastq, .fastq.gz o .fq.gz, con secuencias de longitud corta (hasta 200pb) del tipo paired-end, es decir un archivo forward y un reverse (2x150pb), que corresponden ambos extremos de un fragmento, provenientes de las plataformas de secuenciamiento Illumina, BGI, Ion Torrent y Oxford Nanopore.

**5. Convención (insertar imagen de convención que correspondea la figura 4 de la tesis**

Why this is important?
Incorporar diferentes programas dirigidos a que permitan estudiar el microbioma de una muestra ambiental de cualquier complejidad, y que cumpla con los principios FAIR, además de velocidad, facilidad de uso, notificación de errores, modularidad y libre uso a los usuarios aún sigue siendo un reto.
Con el objetivo de superar estos desafíos, se propone el flujo de trabajo Metabiome, un pipeline bioinformático producto de una exhaustiva revisión de las diferentes metodologías y programas empleados para el análisis de secuencias metagenómicas. Metabiome fue ideado como un pipeline aplicable por módulos, dado que en ocasiones no se requiere que todos los pasos de un flujo de trabajo sean ejecutados, lo que permite omitir pasos incensarios, y evitar gastar memoria de forma innecesaria. El pipeline permite perfilar la comunidad desde un enfoque taxonómico y funcional basado en secuencias ensambladas o no ensambladas, e incluye un paso de extracción del gen 16S que da la posibilidad de realizar análisis comparativos, en relación a la abundancia, entre la composición bacteriana de una muestra analizada desde diferente enfoques. Para facilitar la lectura del texto, consultar el glosario de términos que se utilizan frecuentemente en las diversas secciones.


.. toctree::
   :maxdepth: 2
   :caption: Contents:



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
