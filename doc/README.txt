Qui ci sono articoli che in svariati modi si occupano di semplificazione mesh. Non sono propriamente greedy perché ho letto che le possibilità greedy sono generalmente piuttosto povere, ma possono darci alcune idee su come variare l'algoritmo per poi testarne l'efficienza.

Mesh Optimization e 'New' Approach ricalcano abbastanza quello che abbiamo fatto noi. Il primo è un articolo importante nell'ambiente, ma appunto non so quanto puuò dirci di nuovo.

Comparison of Simplification Algorithms e Feature-preserving... fanno un excursus sui diversi metodi di semplificazione più o meno canonici.

Di particolare interesse c'è Generalized Contraction che introduce un'estensione dell'edge contraction considerando vertici anche non formalmente uniti da un lato. Però potrebbe complicarci più la vita che altro.. specie se si legge nel dettaglio le casistiche presentate nell'articolo.

La TESI (http://amslaurea.unibo.it/3763/1/TESI.pdf) presenta una tecnica di semplificazione che include anche uno step detto di smoothing in cui si usa il laplaciano sulla mesh con un esponente p<1 per concentrare gli elementi nei punti ad alta curvatura.
In generale, l'idea di smoothing - magari facendo una media locale sui lati ad ogni TOT iterate - è qualcosa che potrebbe aiutare forse. 

Entropy-based simplification presenta un modo alternativo alla definizione della parte geometrica.

Simplification Theory fa giusto uno studio sulla qualità della semplificazione su basi di geometria differenziale, trattando il caso che per noi è OnlyGeo.


COME IDEE PERSONALI:
oltre a fare il collasso simultaneo di N edge suff distanti, si potrebbe
- fare il collasso di tutti quelli sotto una certa soglia delta costo (unico davvero greedy trovato in delle slide)
- non usare il punto ottimale di collasso ma solo quello medio (ma non penso che si guadagni molto)
- non usare neanche il punto medio ma SOLO gli estremi. Questo metodo è molto usato e si chiama half edge contraction o decimazione. Ci permetterebbe di risparmiare su diversi controlli penso perchè i nodi rimarrebbero fissi.
 
