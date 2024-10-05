@echo off
g++ -std=c++17 -O2 -DLOCAL -o run b.cpp
g++ -std=c++17 -O2 -DLOCAL -o slow _.cpp
g++ -std=c++17 -O2 -DLOCAL -o gen gen.cpp
for /l %%x in (1, 1, 10000) do (
    gen > in
    run < in > out
    slow < in > out2
    fc out out2 > fc.txt
    if ERRORLEVEL 1 (
        echo ### FAILED INPUT
        type in
        echo ### OUTPUT
        type out
        echo ### CORRECT OUTPUT
        type out2
        pause
    )
)
echo ### Finished testing.
pause
