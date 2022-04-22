import React from "react";

import { ListItem, Stack, Link, Typography, Divider } from "@mui/material";
import { Link as RouterLink } from "react-router-dom";

const SearchCard = (props) => {
  return (
    <ListItem divider
      alignItems="flex-start"
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
        {/* primary */}
        <Stack direction="row" spacing={2} alignItems="center">
          <Link
            sx={{ fontSize: "24px" }}
            component={RouterLink}
            to="about" // TODO change to Team page
            underline="hover"
          >
           <b>{props.title}</b> 
          </Link>
          {props.primary}
        </Stack>
        {/* secondary */}
        <Typography sx={{ color: 'text.secondary' }}>{props.secondary}</Typography>
        {/* tertiary */}
        <Stack direction="row" alignItems="center" spacing={2}>
          {props.tertiary}
        </Stack>
      </Stack>
    </ListItem>
  );
};

export default SearchCard;
