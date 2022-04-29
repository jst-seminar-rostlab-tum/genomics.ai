import React from "react";

import { List } from "@mui/material";

// display the listing of the searched item
const SearchResultList = (props) => {
  const ListItemWrapper = props.listItemWrapper;

  return (
    <List style={{ display: "flex", flexDirection: "column" }}>
      {props.searchedData.map((searchedItem) => (
        <ListItemWrapper key={searchedItem.id} item={searchedItem} />
      ))}
    </List>
  );
};

export default SearchResultList;
