FROM python:3.8.1-alpine
RUN apk update
RUN apk add make automake gcc g++ git

RUN pip install SEIR

CMD SEIR
