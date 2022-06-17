import React from 'react';

import { List } from '@mui/material';

// display the listing of the searched item
function SearchResultList({
  listItemWrapper, searchResult, user, fetchSearchHandler,
}) {
  const ListItemWrapper = listItemWrapper;
  return (
    <List style={{
      display: 'flex',
      flexDirection: 'column',
    }}
    >
      {searchResult.map((searchedItem) => (
        <ListItemWrapper
          key={searchedItem.id}
          item={searchedItem}
          user={user}
          onAction={fetchSearchHandler}
        />
      ))}
    </List>
  );
}

export default SearchResultList;
