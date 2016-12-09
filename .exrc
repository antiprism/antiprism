if &cp | set nocp | endif
map  :py3f ~/bin/clang-format.py
let s:cpo_save=&cpo
set cpo&vim
vmap gx <Plug>NetrwBrowseXVis
nmap gx <Plug>NetrwBrowseX
vnoremap <silent> <Plug>NetrwBrowseXVis :call netrw#BrowseXVis()
nnoremap <silent> <Plug>NetrwBrowseX :call netrw#BrowseX(expand((exists("g:netrw_gx")? g:netrw_gx : '<cfile>')),netrw#CheckIfRemote())
inoremap  
imap  :py3f ~/bin/clang-format.pyi
inoremap # X#
cabbr nall :n *.c *.cc *.cpp *.h 
let &cpo=s:cpo_save
unlet s:cpo_save
set autoindent
set autowrite
set backspace=2
set backup
set backupdir=~/.vim
set cindent
set completeopt=menuone,longest,preview
set expandtab
set fileencodings=ucs-bom,utf-8,default,latin1
set helplang=en
set hidden
set laststatus=2
set listchars=eol:-
set pastetoggle=<F12>
set printoptions=paper:a4
set ruler
set runtimepath=~/.vim,~/.vim/bundle/Vundle.vim,~/.vim/bundle/syntastic,/var/lib/vim/addons,/usr/share/vim/vimfiles,/usr/share/vim/vim74,/usr/share/vim/vimfiles/after,/var/lib/vim/addons/after,~/.vim/after,~/.vim/bundle/Vundle.vim,~/.vim/bundle/Vundle.vim/after,~/.vim/bundle/syntastic/after
set shiftwidth=3
set softtabstop=2
set suffixes=.bak,~,.swp,.o,.info,.aux,.log,.dvi,.bbl,.blg,.brf,.cb,.ind,.idx,.ilg,.inx,.out,.toc
set tags=./tags,tags
" vim: set ft=vim :
