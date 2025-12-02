# example usage:


    docker run --rm -it -v `pwd`:/data us-central1-docker.pkg.dev/methods-dev-lab/lraa/cas \
         --matrix-dir /data/refQuant^gene-sparseM-forCAS \
         --api-token ${PRIVATE_TOKEN} \
         --output-prefix /data/cas.test.out


    