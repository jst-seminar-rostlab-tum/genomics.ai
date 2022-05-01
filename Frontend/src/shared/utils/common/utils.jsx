import { Redirect } from 'react-router-dom';
import React from 'react';

export default function removeItemFromArray(arr, value) {
  const index = arr.indexOf(value);
  if (index > -1) {
    arr.splice(index, 1);
  }
  return arr;
}

export function guardedPage(page) {
  if (!!localStorage.jwt && !!localStorage.user) {
    return page;
  }
  return <Redirect to="/" />;
}

export function getAuthHeader() {
  return {
    auth: localStorage.getItem('jwt'),
  };
}

export function getAuthAndJsonHeader() {
  return {
    auth: localStorage.getItem('jwt'),
    'content-type': 'application/json',
  };
}

export function createUrl(path, pathParam, value) {
  return path.replace(`:${pathParam}`, value);
}

export function setSeachCategoryInUrl(path, value) {
  return createUrl(path, 'searchCategory', value);
}
