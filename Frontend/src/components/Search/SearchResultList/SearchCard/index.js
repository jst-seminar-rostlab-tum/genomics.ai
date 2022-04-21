import React from "react";

import { ListItem, Stack } from "@mui/material";

const SearchCard = (props) => {
  return (
    <ListItem
      alignItems="flex-start"
      key={props.key}
      secondaryAction={props.action}
    >
      {
        // TODO add avatar if exists
        /* <ListItemAvatar>
        <Avatar
          sx={{ bgcolor: randomColor(), width: 45, height: 45 }}
          alt="Remy Sharp"
          src="/static/images/avatar/1.jpg"
        />
      </ListItemAvatar> */
      }
      <Stack direction="column" spacing={0.5}>
        <Stack direction="row" spacing={2}>
          {props.primary}
        </Stack>
        <div>{props.secondary}</div>
        <Stack direction="row" alignItems="center" spacing={2}>
          {props.tertiary}
        </Stack>
      </Stack>
    </ListItem>
  );
};

export default SearchCard;
