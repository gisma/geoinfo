# ---- Berechnung change from - to ----
# Die Darstellung der Veränderungen von jeder Klasse zu jeder Klasse ist mit
# "Bordmitteln"  etwas aufwendiger. Da jede Kategorie mit jeder verglichen
# werden soll müssen dazu ein paar Hilfsfunktionen und Schleifen (loops)
# eingesetzt werden Hilfsfunktion changefrom() vergleicht je zwei Raster auf die
# Kategorien i,j
changefrom=function(r1,r2,i,j){r1==i & r2==j}
