FROM python:3.5-slim

WORKDIR /data

RUN pip3 install matplotlib \
                 scipy

COPY *.py /usr/local/bin/
COPY Atha.tsv /data 

CMD bash
