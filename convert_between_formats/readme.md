## Visualization of cool files

Use HiGlass to visualize cool files. The simple example below stops any running higlass containers, removes them, pulls the latest version and runs it:

```bash
    # Do once
    docker pull higlass/higlass-docker:v0.6.1 # higher versions are experimental and may or may not work

    docker stop higlass-container;
    docker rm higlass-container;

    docker run --detach \
               --publish 8989:80 \
               --volume /Users/sreenivasan/Documents/Works/janaLmnb1_public_data/scHiC/test:/data \
               --volume /tmp:/tmp \
               --name higlass-container \
             higlass/higlass-docker:v0.6.1
```

The higlass website should now be visible at http://localhost:8989. Take a look at the documentation for adding a new track to see how to display data.

To load a cooler file, first convert to mcool to visualize using HiGlass

```bash
# do once
mamba create --name cooler python=3.12
conda activate cooler
pip install cooler

pip install higlass-manage

# Convert the cool file to mcool file
cooler zoomify --balance DNMT3Awtc13.allValidPairs.hic_10000_10000.cool

docker exec higlass-container python higlass-server/manage.py \
   ingest_tileset \
   --filename /data/DNMT3Awtc13.allValidPairs.hic_10000_10000.mcool \
   --datatype matrix \
   --filetype cooler
```

Note: HiPiler build on higlass to make visualization and book marking even better
