FROM python:3.6

RUN pip install matplotlib numpy pandas seaborn statsmodels scipy

RUN git clone https://github.com/tuvangezer/selfish/