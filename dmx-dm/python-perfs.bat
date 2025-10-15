echo Set Python environment variables...

@echo off
call "%ProgramFiles%\GNU\WPy-3662\scripts\env.bat

python.exe noiseAnalysis-e2e.py

pause
