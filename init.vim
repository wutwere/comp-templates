call plug#begin()

Plug 'itchyny/lightline.vim'
let g:lightline = { 'colorscheme': 'github' }

Plug 'sheerun/vim-polyglot'

Plug 'cormacrelf/vim-colors-github'
Plug 'tomasr/molokai'
Plug 'dracula/vim', { 'as': 'dracula' }

Plug 'ThePrimeagen/vim-be-good'

call plug#end()

autocmd VimEnter * if len(filter(values(g:plugs), '!isdirectory(v:val.dir)'))
  \| PlugInstall --sync | source $MYVIMRC
\| endif

cd ~/Desktop/code

syntax on
color github

set background=light
set termguicolors
set nu
set rnu
set expandtab
set tabstop=2
set shiftwidth=2
set mouse=a

inoremap {<CR> {<CR>}<Esc>O

let g:template =<< END
#include <bits/stdc++.h>
using namespace std;

#ifdef LOCAL
#include "templates/debug.h"
#else
#define dbg(x...)
#endif

// make it simpler

int main() {
  cin.tie(0)->sync_with_stdio(0);
  for (int __, _ = (cin >> __, 0); ++_ <= __;) {
    cout << "Case #" << _ << ": ";
  }
}
END

function Run_clipboard()
  call writefile(split(getreg('+'), '\n'), 'in')
  if expand('%:e') == 'cpp'
    execute '!start cmd /c "run < in & pause"'
  elseif expand('%:e') == 'py'
    execute 'w | !start cmd /c "py ' . expand('%:t') . ' < in & pause"'
  endif
endfunction

autocmd filetype cpp nnoremap <C-N> :<C-U>%d \| call setline(1, g:template)<CR>G2k
autocmd filetype cpp nnoremap <C-B> :<C-U>w \| !g++ -O2 -DLOCAL -std=c++17 %:r.cpp -o run<CR>
nnoremap <C-C> :<C-U>call Run_clipboard()<CR><CR>
autocmd filetype cpp nnoremap <C-S> :<C-U>!start cmd /c "run & timeout /t -1 /nobreak"<CR><CR>
nnoremap <C-A> :<C-U>%y+<CR>
