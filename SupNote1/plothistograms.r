G4Plus2kflanks = read.table('GQuadPlus2kflanks_chr21Ends.tab', header = TRUE)
Left = G4Plus2kflanks[which(G4Plus2kflanks$event_relative_pos=='L'),]
Right = G4Plus2kflanks[which(G4Plus2kflanks$event_relative_pos=='R'),]
Center = G4Plus2kflanks[which(G4Plus2kflanks$event_relative_pos=='C'),]
Center$event_coord_inwindow = round(Center$event_coord_inwindow / (Center$motif_end - Center$motif_start +1), digits = 3)

LeftEnds <- as.numeric(unlist(Left$event_coord_inwindow))
png(filename="LeftEnds.png")
hist(-LeftEnds, breaks = 2000, ylim = c(0,200))
dev.off

CenterEnds <- as.numeric(unlist(Center$event_coord_inwindow))
png(filename="CenterEnds.png")
hist(CenterEnds, breaks = 100, ylim = c(0,200))
dev.off

RightEnds <- as.numeric(unlist(Right$event_coord_inwindow))
png(filename="RightEnds.png")
hist(RightEnds, breaks = 2000, ylim = c(0,200))
dev.off

G4Plus2kflanks = read.table('GQuadPlus2kflanks2_chr21Ends.tab', header = TRUE)
Left = G4Plus2kflanks[which(G4Plus2kflanks$event_relative_pos=='L'),]
Right = G4Plus2kflanks[which(G4Plus2kflanks$event_relative_pos=='R'),]
Center = G4Plus2kflanks[which(G4Plus2kflanks$event_relative_pos=='C'),]
Center$event_coord_inwindow = round(Center$event_coord_inwindow / (Center$motif_end - Center$motif_start +1), digits = 3)

LeftEnds <- as.numeric(unlist(Left$basephred))
hist(-LeftEnds, breaks = 2000, ylim = c(0,200))

CenterEnds <- as.numeric(unlist(Center$basephred))
hist(CenterEnds, breaks = 100, ylim = c(0,200))

RightEnds <- as.numeric(unlist(Right$basephred))
hist(RightEnds, breaks = 2000, ylim = c(0,200))