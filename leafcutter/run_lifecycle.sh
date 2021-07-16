#!/usr/bin/env bash

LOCAL_NAME=leafcutter
DIR_NAME=leafcutter

RKSTR8_REL_PATH=../../rkstr8
RKSTR8_REQS_REL_PATH=../../requirements.small.txt

MASTER_SCRIPT=run_leafcutter.py

AWS_CONFIG_LOCAL=config

USER_NAME=$1
REPOSITORY_NAME=leafcutter
REPOSITORY_TAG=latest

REMOTE_PATH="${USER_NAME}/${REPOSITORY_NAME}:${REPOSITORY_TAG}"

# Make config

echo -e "[default]\nregion=us-east-1\n" > $AWS_CONFIG_LOCAL

if [ -e "$AWS_CONFIG_LOCAL" ]; then
    echo "config exists!"
fi

cp -R $RKSTR8_REL_PATH .
cp $RKSTR8_REQS_REL_PATH requirements.txt

cp $MASTER_SCRIPT .

echo "Building..."
docker build --tag=${LOCAL_NAME} .
echo "Tagging..."
docker tag ${LOCAL_NAME} ${REMOTE_PATH}
echo "Pushing..."
docker push ${REMOTE_PATH}

# Tear down rkstr8 and requirements
my_dir=$(basename `pwd`)
if [ "$my_dir" = "$DIR_NAME" ]; then
    rm -rf rkstr8
    rm requirements.txt
    rm config
fi
