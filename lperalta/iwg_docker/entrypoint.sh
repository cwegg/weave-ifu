#!/bin/bash

if [ -n "$HOST_USER_ID" ] && [ -n "$HOST_USER_NAME" ] && [ -n "$HOST_GROUP_ID" ] && [ -n "$HOST_GROUP_NAME" ]
then
    echo "Creating user $HOST_USER_NAME:$HOST_GROUP_NAME ($HOST_USER_ID:$HOST_GROUP_ID)"
    
    groupadd -g ${HOST_GROUP_ID} ${HOST_GROUP_NAME}
    useradd -m -s /bin/bash -u ${HOST_USER_ID} -g ${HOST_GROUP_ID} ${HOST_USER_NAME}
    
    chown ${HOST_USER_NAME}:${HOST_GROUP_NAME} /data /data/*
else
    echo "Skiping user creation because the following environment variables are not defined: HOST_USER_ID HOST_USER_NAME HOST_GROUP_ID HOST_GROUP_NAME"
fi

exec "$@"

