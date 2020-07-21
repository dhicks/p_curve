PAPER := paper
SCRIPT := scripts

all: script paper

paper: 
	$(MAKE) -C $(PAPER)

script:
	$(MAKE) -C $(SCRIPT)

