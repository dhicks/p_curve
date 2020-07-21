PAPER := paper
SCRIPT := scripts

all: script ms

script:
	$(MAKE) -C $(SCRIPT)

ms:
	$(MAKE) -C $(PAPER)
