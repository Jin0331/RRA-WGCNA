library(RCy3)
library(igraph)

.defaultBaseUrl <- "http://192.168.0.7:1234/v1"
cytoscapePing(base.url = .defaultBaseUrl)
cytoscapeVersionInfo(base.url = .defaultBaseUrl)


sif <- system.file("extdata","galFiltered.sif",package="RCy3")
importNetworkFromFile(sif, base.url = .defaultBaseUrl)

csv <- system.file("extdata","galExpData.csv", package="RCy3")
data <- read.csv(csv, stringsAsFactors = FALSE)

loadTableData(data, data.key.column="name", base.url = .defaultBaseUrl)
createNetworkFromIgraph(ig,"myIgraph", base.url = .defaultBaseUrl)

ig2 <- createIgraphFromNetwork("myIgraph", base.url = .defaultBaseUrl)


getNetworkList(, base.url = .defaultBaseUrl)
