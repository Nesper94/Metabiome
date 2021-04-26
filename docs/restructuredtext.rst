Sphinx y reSTructuredtext
===========================

*Sphinx* es un generador de documentación, que utiliza el lenguaje **restructuredtext** (una sola palabra).

* Lenguaje de marcado simple y discreto.

sphinx-quickstart
---
`Guía rápida de instalación y características del contenido`__ <https://sphinx-tutorial.readthedocs.io/start/>__

_________________________________________________________________

Sintaxis básica
===============


* **Títulos o cabeceras de título **:Se indican empleando un signo de puntuación debajo de la cabecera, esto me permite estructuración del texto:

Título 1
========

Título 2
-------

Título 3
^^^^^^^^

Título 4
""""""""

* Algunos formatos: **negrita**, *cursiva*,
``código``.

    ``Todo lo que quiera que quede monoespaciado va dentro de doble comilla invertida``

* Barra invertida (\\) para anular el significado especial. Para obtener una barra invertida literal, usar dible barra invertida ("\\"). 

*ejemplo* **barra**  ``invertida``

\*ejemplo* \**barra**  \``invertida``



* **Hipervínculos** 
    + Hipervínculos con destinos explicitos:

`Git <https://git-scm.com>`_


https://www.python.org/


`Sphinx`__

__ http://www.sphinx-doc.org



* **Texto preformateado**
¿Dónde se separa el texto y cómo se rompen las líneas?


Por ejemplo, si se supone que puse esa pregunta en otra línea es poque quería que se escribiera en una línea separadas, pero no pasó así, ¿por qué?

1. HTML no reconoce el espacio que se agrega al código, incluida la barra espaciadora. Si coloca veinte espacios entre una palabra y la palabra que le sigue, el navegador mostrará un solo espacio allí. Esto significa "**Colapso de espacios en blanco**".
2. Para tener más control sobre cómo se lee el documento debo agregar texto **preformateado**.
3. Dos líneas de separación indican un nuevo párrafo.

Ahora se supone que es un nuevo párrafo porqué tiene dos lineas de espacio

* **Literal blocks o "bloques literales":** 
    Muestran exactamente lo que escribo en texto plano:

:: 

    Un párrafo que contiene solo dos doble punto indica que el siguiente texto indentado o      
    entre comillas es un bloque literal. *como esto* que debería quedar en cursiva, o
    **esto**, que debería mostrarse en negrita. 
    
    
Este es ejemplo del uso de dos puntos dobles seguido de texto:: 

    *todo lo que* está indentado se muestra literal, si al final del párrafo anterior hay un "::" sin espacio. 

* **Line blocks o "bloque de línea"
    
    Dar forma de versos o adornar el texto

| Este es un ejemplo de bloques de línea
| pero parece que no funciona
| cada línea empiea con una barra vertical





Table of Contents Tree
----------------------


**toctree**
    Sphinx permite combinar varias páginas en una jerarquía cohesiva. La directiva toctree es parte fundamental de esta estructura. 
 ____________________________________________

* **Ejemplos de código:** Para mostrar ejemplos de código, o cualquier texto que no deba formatearse, finalice el párrafo antes del ejemplo de código con dos dos puntos.::

``This is the paragraph preceding the code sample::

    Se supone que todo lo que yo escriba acá se debe ver como código

















