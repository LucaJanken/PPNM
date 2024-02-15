#ifndef HAVE_GENLIST_H
#define HAVE_GENLIST_H

#include <iostream>
#include <cassert>

template <typename T>
class genlist {
private:
    struct Node {
        T item;
        Node* next;
        Node(const T& item, Node* next = nullptr) : item(item), next(next) {}
    };

    Node* head;
    Node* tail;
    size_t _size;

public:
    const size_t& size = _size;

    genlist() : head(nullptr), tail(nullptr), _size(0) {}

    ~genlist() {
        Node* current = head;
        while (current != nullptr) {
            Node* next = current->next;
            delete current;
            current = next;
        }
    }

    genlist(const genlist& other) : head(nullptr), tail(nullptr), _size(0) {
        Node* current = other.head;
        while (current != nullptr) {
            add(current->item);
            current = current->next;
        }
    }

    genlist& operator=(const genlist& other) {
        if (this != &other) {
            while (head != nullptr) {
                Node* toDelete = head;
                head = head->next;
                delete toDelete;
            }
            head = tail = nullptr;
            _size = 0;

            Node* current = other.head;
            while (current != nullptr) {
                add(current->item);
                current = current->next;
            }
        }
        return *this;
    }

    T& operator[](size_t i) {
        assert(i < _size);
        Node* current = head;
        for (size_t index = 0; index < i; ++index) {
            current = current->next;
        }
        return current->item;
    }

    void add(const T& item) {
        Node* newNode = new Node(item);
        if (head == nullptr) {
            head = tail = newNode;
        } else {
            tail->next = newNode;
            tail = newNode;
        }
        _size++;
    }
};

#endif // HAVE_GENLIST_H
