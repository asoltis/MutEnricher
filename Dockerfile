FROM python:3.7.7

ARG BUILD_DATE

LABEL maintainer="anthonyrsoltis@gmail.com"
LABEL version="1.3.3"
LABEL build-date=$BUILD_DATE
LABEL description="MutEnricher: somatic coding and noncoding mutation enrichment analysis tool."
LABEL url="https://github.com/asoltis/MutEnricher"

WORKDIR /usr/src/app

RUN curl -L http://github.com/asoltis/mutenricher/archive/master.tar.gz -o MutEnricher-master.tar.gz

RUN tar -xzvf MutEnricher-master.tar.gz

RUN rm -f MutEnricher-master.tar.gz

RUN mkdir MutEnricher

RUN cp -r MutEnricher-master/* MutEnricher/

RUN pip install numpy==1.17.3

RUN pip install scipy==1.3.1

RUN pip install cython==0.29.13

RUN pip install pysam==0.15.2

RUN pip install cyvcf2==0.20.0

WORKDIR /usr/src/app/MutEnricher/math_funcs

RUN python setup.py build_ext --inplace

WORKDIR /usr/src/app/MutEnricher

CMD ["python", "mutEnricher.py"]

