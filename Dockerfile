FROM 812206152185.dkr.ecr.us-west-2.amazonaws.com/latch-base:fe0b-main

COPY requirements.txt /root/requirements.txt
RUN pip install --no-cache-dir -r /root/requirements.txt

COPY wf /root/wf

ARG tag
ENV FLYTE_INTERNAL_IMAGE $tag
WORKDIR /root
