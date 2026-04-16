# zebrafish-ros-normalization-online

Aplicación web estática para normalizar y analizar datos de fluorescencia DCF en embriones de pez cebra directamente en el navegador usando Pyodide.

## Qué hace

- corre completamente del lado del cliente;
- no requiere servidor ni backend;
- acepta múltiples archivos `.csv`;
- ignora automáticamente archivos que no sean `.csv`;
- procesa cada archivo por separado;
- normaliza por fecha usando el DMSO del mismo día;
- detecta outliers por `fecha × condición`;
- produce ramas paralelas `with_outliers` y `without_outliers`;
- muestra gráficas interactivas y warnings;
- permite descargar un ZIP con tablas de salida por archivo.

## Estructura

- `index.html`: interfaz principal
- `style.css`: estilo visual, inspirado en `Clonalyzer-2`
- `app.js`: lógica de interfaz, Pyodide, Plotly y descargas
- `zebrafish_ros_engine.py`: motor analítico en Python para navegador

## Desarrollo local

Para vista rápida local:

```bash
python3 -m http.server 8000
```

Luego abre:

```text
http://localhost:8000
```

## GitHub Pages

El proyecto está diseñado para desplegarse como sitio estático en GitHub Pages.

Pasos típicos:

1. publicar `main` en GitHub;
2. activar GitHub Pages desde la rama `main` y carpeta raíz;
3. esperar a que GitHub publique el sitio.

## Notas

- La primera carga de Pyodide puede tardar varios segundos.
- La app descarga dependencias desde CDN.
- El procesamiento ocurre completamente en el navegador del usuario.
