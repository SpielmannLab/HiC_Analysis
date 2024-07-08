## Visualization of cool files
Use HiGlass to visualize cool files. The simple example below stops any running higlass containers, removes them, pulls the latest version and runs it:

    # Do once
    docker pull higlass/higlass-docker:v0.6.1 # higher versions are experimental and may or may not work

    docker stop higlass-container;
    docker rm higlass-container;

    docker run --detach \
               --publish 8989:80 \
               --volume ~/Desktop/2024.05_cool_format/raw_cool:/data \
               --volume ~/tmp:/tmp \
               --name higlass-container \
             higlass/higlass-docker:v0.6.1

The higlass website should now be visible at http://localhost:8989. Take a look at the documentation for adding a new track to see how to display data.

To load a cooler file, first convert to mcool to visualize using HiGlass

    # do once
    mamba create --name cooler python=3.7
    conda activate cooler
    pip install cooler

    # Convert the cool file to mcool file
    cooler zoomify --balance DNMT3Awtc13.allValidPairs.hic_10000_10000.cool

    docker exec higlass-container python higlass-server/manage.py \
      ingest_tileset \
      --filename /data/DNMT3Awtc13.allValidPairs.hic_10000_10000.mcool \
      --datatype matrix \
      --filetype cooler
