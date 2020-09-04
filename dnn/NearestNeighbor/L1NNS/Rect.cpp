#include"Rect.h"

Rect::Rect() {}

Rect::Rect(double x1, double x2, double y1, double y2)
{
	this->l = x1;
	this->r = x2;
	this->d = y1;
	this->u = y2;
	this->ld = new Point(x1, y1);
	this->lu = new Point(x1, y2);
	this->rd = new Point(x2, y1);
	this->ru = new Point(x2, y2);
}

Rect::~Rect() {
	delete this->ld;
	delete this->lu;
	delete this->rd;
	delete this->ru;
}

double Rect::getl() {
	return this->l;
}

double Rect::getr() {
	return this->r;
}

double Rect::getu() {
	return this->u;
}

double Rect::getd() {
	return this->d;
}

Point* Rect::getld() {
	return this->ld;
}

Point* Rect::getlu() {
	return this->lu;
}

Point* Rect::getrd() {
	return this->rd;
}

Point* Rect::getru() {
	return this->ru;
}

int Rect::getid() {
	return this->id;
}

void Rect::setid(int _id) {
	this->id = _id;
}