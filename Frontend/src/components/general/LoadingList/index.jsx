import React from 'react';
import CircularProgress from '@mui/material/CircularProgress';
import SpaceBetweenCards from 'components/general/SpaceBetweenCards';

/**
 * @param {isLoading} whether the data is still loading (spinner will be shown)
 * @param {elements} a list of list elements to iterate over
 *                   - each element needs to have the id property
 * @param {cardBuilder} takes an element and returns a card
 * @param {noElementsMessage} a message to show if there are no elements
 * @returns react elements
 */
function LoadingList({
  isLoading, elements, cardBuilder, noElementsMessage,
}) {
  function loadingWrapper(orElse) {
    return isLoading ? (
      <CircularProgress />
    ) : orElse;
  }

  return (
    loadingWrapper(
      elements.length ? (
        elements.map((elem) => (
          <div key={elem.id}>
            {cardBuilder(elem)}
            <SpaceBetweenCards />
          </div>
        ))
      ) : (
        <span>{noElementsMessage}</span>
      ),
    )
  );
}

export default LoadingList;
