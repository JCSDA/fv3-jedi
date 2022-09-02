## How to generate fv3grid files for other resolutions?

- Compile [fv3-bundle](https://github.com/JCSDA/fv3-bundle)
- Update **npx** and **npy** in the **geometry_gfs** test:

    cd /path/to/your/build/fv3-bundle/fv3-jedi
    vi test/testinput/geometry_gfs.yaml

- Run the **geometry_gfs** test:

    ctest -R fv3jedi_test_tier1_geometry_gfs

- fv3grid files are located in the **test** directory
