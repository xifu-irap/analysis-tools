echo Set Python environment variables...

@echo off
call C:\GNU\WPy64-313110\scripts\env.bat

"C:\GNU\WPy64-313110\python\python.exe"  fdbkDelayAnalysis_e2e.py

pause
