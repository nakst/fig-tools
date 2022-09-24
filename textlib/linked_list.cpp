// This file is part of the Essence operating system.
// It is released under the terms of the MIT license -- see LICENSE.md.
// Written by: nakst.

template <class T>
struct LinkedList;

template <class T>
struct LinkedItem {
	void RemoveFromList();

	LinkedItem<T> *previousItem;
	LinkedItem<T> *nextItem;

	// TODO Separate these out?
	struct LinkedList<T> *list;
	T *thisItem;
};

template <class T>
struct LinkedList {
	void InsertStart(LinkedItem<T> *item);
	void InsertEnd(LinkedItem<T> *item);
	void InsertBefore(LinkedItem<T> *newItem, LinkedItem<T> *beforeItem);
	void Remove(LinkedItem<T> *item);

	void Validate(int from); 

	LinkedItem<T> *firstItem;
	LinkedItem<T> *lastItem;

	size_t count;

#ifdef DEBUG_BUILD
	bool modCheck;
#endif
};

struct SimpleList {
	void Insert(SimpleList *link, bool start);
	void Remove();

	union { SimpleList *previous, *last; };
	union { SimpleList *next, *first; };
};

#ifndef SHARED_DEFINITIONS_ONLY

template <class T>
void LinkedItem<T>::RemoveFromList() {
	if (!list) {
		EsPanic("LinkedItem::RemoveFromList - Item not in list.\n");
	}

	list->Remove(this);
}

template <class T>
void LinkedList<T>::InsertStart(LinkedItem<T> *item) {
#ifdef DEBUG_BUILD
	if (modCheck) EsPanic("LinkedList::InsertStart - Concurrent modification\n");
	modCheck = true; EsDefer({modCheck = false;});
#endif

	if (item->list == this) EsPanic("LinkedList::InsertStart - Inserting an item that is already in this list\n");
	if (item->list) EsPanic("LinkedList::InsertStart - Inserting an item that is already in a list\n");

	if (firstItem) {
		item->nextItem = firstItem;
		item->previousItem = nullptr;
		firstItem->previousItem = item;
		firstItem = item;
	} else {
		firstItem = lastItem = item;
		item->previousItem = item->nextItem = nullptr;
	}

	count++;
	item->list = this;
	Validate(0);
}

template <class T>
void LinkedList<T>::InsertEnd(LinkedItem<T> *item) {
#ifdef DEBUG_BUILD
	if (modCheck) EsPanic("LinkedList::InsertEnd - Concurrent modification\n");
	modCheck = true; EsDefer({modCheck = false;});
#endif

	if (item->list == this) EsPanic("LinkedList::InsertEnd - Inserting a item that is already in this list\n");
	if (item->list) EsPanic("LinkedList::InsertEnd - Inserting a item that is already in a list\n");

	if (lastItem) {
		item->previousItem = lastItem;
		item->nextItem = nullptr;
		lastItem->nextItem = item;
		lastItem = item;
	} else {
		firstItem = lastItem = item;
		item->previousItem = item->nextItem = nullptr;
	}

	count++;
	item->list = this;
	Validate(1);
}

template <class T>
void LinkedList<T>::InsertBefore(LinkedItem<T> *item, LinkedItem<T> *before) {
#ifdef DEBUG_BUILD
	if (modCheck) EsPanic("LinkedList::InsertBefore - Concurrent modification\n");
	modCheck = true; EsDefer({modCheck = false;});
#endif

	if (item->list == this) EsPanic("LinkedList::InsertBefore - Inserting a item that is already in this list\n");
	if (item->list) EsPanic("LinkedList::InsertBefore - Inserting a item that is already in a list\n");

	if (before != firstItem) {
		item->previousItem = before->previousItem;
		item->previousItem->nextItem = item;
	} else {
		firstItem = item;
		item->previousItem = nullptr;
	}

	item->nextItem = before;
	before->previousItem = item;

	count++;
	item->list = this;
	Validate(3);
}

template <class T>
void LinkedList<T>::Remove(LinkedItem<T> *item) {
#ifdef DEBUG_BUILD
	if (modCheck) EsPanic("LinkedList::Remove - Concurrent modification\n");
	modCheck = true; EsDefer({modCheck = false;});
#endif

	if (!item->list) EsPanic("LinkedList::Remove - Removing an item that has already been removed\n");
	if (item->list != this) EsPanic("LinkedList::Remove - Removing an item from a different list (list = %x, this = %x)\n", item->list, this);

	if (item->previousItem) {
		item->previousItem->nextItem = item->nextItem;
	} else {
		firstItem = item->nextItem;
	}

	if (item->nextItem) {
		item->nextItem->previousItem = item->previousItem;
	} else {
		lastItem = item->previousItem;
	}

	item->previousItem = item->nextItem = nullptr;
	item->list = nullptr;
	count--;
	Validate(2);
}

template <class T>
void LinkedList<T>::Validate(int from) {
#ifdef DEBUG_BUILD
	if (count == 0) {
		if (firstItem || lastItem) {
			EsPanic("LinkedList::Validate (%d) - Invalid list (1)\n", from);
		}
	} else if (count == 1) {
		if (firstItem != lastItem
				|| firstItem->previousItem
				|| firstItem->nextItem
				|| firstItem->list != this
				|| !firstItem->thisItem) {
			EsPanic("LinkedList::Validate (%d) - Invalid list (2)\n", from);
		}
	} else {
		if (firstItem == lastItem
				|| firstItem->previousItem
				|| lastItem->nextItem) {
			EsPanic("LinkedList::Validate (%d) - Invalid list (3) %x %x %x %x\n", from, firstItem, lastItem, firstItem->previousItem, lastItem->nextItem);
		}

		{
			LinkedItem<T> *item = firstItem;
			size_t index = count;

			while (--index) {
				if (item->nextItem == item || item->list != this || !item->thisItem) {
					EsPanic("LinkedList::Validate (%d) - Invalid list (4)\n", from);
				}

				item = item->nextItem;
			}

			if (item != lastItem) {
				EsPanic("LinkedList::Validate (%d) - Invalid list (5)\n", from);
			}
		}

		{
			LinkedItem<T> *item = lastItem;
			size_t index = count;

			while (--index) {
				if (item->previousItem == item) {
					EsPanic("LinkedList::Validate (%d) - Invalid list (6)\n", from);
				}

				item = item->previousItem;
			}

			if (item != firstItem) {
				EsPanic("LinkedList::Validate (%d) - Invalid list (7)\n", from);
			}
		}
	}
#else
	(void) from;
#endif
}

void SimpleList::Insert(SimpleList *item, bool start) {
	if (item->previous || item->next) {
		EsPanic("SimpleList::Insert - Bad links in %x.\n", this);
	}

	if (!first && !last) {
		item->previous = this;
		item->next = this;
		first = item;
		last = item;
	} else if (start) {
		item->previous = this;
		item->next = first;
		first->previous = item;
		first = item;
	} else {
		item->previous = last;
		item->next = this;
		last->next = item;
		last = item;
	}
}

void SimpleList::Remove() {
	if (previous->next != this || next->previous != this) {
		EsPanic("SimpleList::Remove - Bad links in %x.\n", this);
	}

	if (previous == next) {
		next->first = nullptr;
		next->last = nullptr;
	} else {
		previous->next = next;
		next->previous = previous;
	}

	previous = next = nullptr;
}

#endif
