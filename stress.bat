@echo off
g++ -std=c++17 -O2 -DLOCAL -o run f.cpp
g++ -std=c++17 -O2 -DLOCAL -o slow _.cpp
g++ -std=c++17 -O2 -DLOCAL -o gen gen.cpp
for /l %%x in (1, 1, 10000) do (
    gen > in.txt
    run < in.txt > out.txt
    slow < in.txt > out2.txt
    fc out.txt out2.txt > fc.txt
    if ERRORLEVEL 1 (
        echo ### FAILED INPUT
        type in.txt
        echo ### OUTPUT
        type out.txt
        echo ### CORRECT OUTPUT
        type out2.txt
        pause
    )
)
echo ### Finished testing.
pause
