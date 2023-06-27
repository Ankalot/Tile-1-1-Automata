# Tile(1, 1) automata
![2023-06-25 22-25-59 (online-video-cutter com)](https://github.com/Ankalot/Tile-1-1-Automaton/assets/74981478/b93f7f76-de85-48ee-b716-9016e3791ba2)  
Cellular automata, but instead of cells there are tiles Tile(1, 1).  
I took theory and code from this repo: https://github.com/shrx/spectre, and modified it for cellular automata.  

# How to use
1) Run `./compile.bat` (or execute commands manually from it).  
You can change `N_ITERATIONS` in `generateTiles.cpp` for different amounts of tiles.  
Information about tiles will be generated in `cells.txt`.  
2) Open `index.html` in your browser.  
3) Click `Choose file` button and select `cells.txt`.  
4) Click the button `Init`.  
5) Now you can click on tiles, marking them `alive`.  
6) Click the button `Start automata`.  

You can click on tiles while the automata is working.  
For a reset, you need to update the page.  
For changing automata rules, you need to change `start` function in `index.html`. 
`numAlive(tiles[i]["neighbours"])` means number of alive neighboring tiles, `numAlive(tiles[i]["vertNeighbours"])` 
means number of alive neighboring tiles that have only a common vertex with the given one.  
