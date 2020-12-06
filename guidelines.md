# Guidelines

- All packages should be installed within a Conda environment, so no
script should install anything unless necessary.

- When development work is done, Conda environments must be exported to a
 **YAML** file using `conda env export --from-history --name env-name > env-name.yml` so that
they can be recreated on installation. It would also be useful to have a
file with platform-dependent information specifying the version of each
package. For example with an environment named `humann2`:

```
conda env export --from-history --name humann2 > version-humann2.yml && echo "# Platform $(lsb_release -d)" >> version-humann2.yml
```

- Remember the [KISS principle](https://en.wikipedia.org/wiki/KISS_principle)
and make things as simple as possible. As Nobel Prize for physics Paul Dirac
 once stated:
> The aim of science is to make difficult things understandable in a simpler
 way.

- A script should stop execution if called with no arguments:
```bash
if [[ "$#" == 0 ]]; then
    echo "No arguments given."
    usage
    exit 1
fi
```

- They also should return some useful info:
```bash
# Useful output info
echo "Input directory: ${input_dir:?'Input directory not set'}"
echo "Output directory: ${out_dir:?'Output directory not set'}"
echo "Number of threads: ${threads:=4}"
echo "<Package> version: $(<package> --version)"
```

- And they should perform basic input validation:
```bash
# Verify that input directory exists
if [ ! -d "$input_dir" ]; then
   echo "$0: Error: $input_dir is not a valid directory."
   exit 1
fi
# Create output directory if it doesn't exists.
if [[ ! -d "$out_dir" ]]; then
    mkdir "$out_dir"
fi
```

Feel free to copy and paste the code in this file when creating your scripts.
