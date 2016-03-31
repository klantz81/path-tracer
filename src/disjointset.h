#include <vector>

struct node {
	node *parent;
	int i, rank;
};
	
class cDisjointSet {
  private:
	std::vector<node *> nodes;
	int elements, sets;

  protected:
    
  public:
	cDisjointSet();
	~cDisjointSet();

	node* MakeSet(int i);
	node* Find(node* a);
	void Union(node* a0, node* a1);

	int ElementCount();
	int SetCount();

	int Reduce();
	void Reset();
};