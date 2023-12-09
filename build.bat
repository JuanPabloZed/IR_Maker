@echo off
cls 

if not exist "venv" (
    python -m venv venv
    ./venv/Scripts/activate.bat
    
    echo Virtual environment created
    echo Installing python libraries

    pip install -r requirements.txt

    echo All libraries installed

)

pyinstaller --noconfirm --onefile --windowed ./IRMaker.py