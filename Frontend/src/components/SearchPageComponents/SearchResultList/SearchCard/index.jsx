import React from 'react';

import {
  ListItem,
  Stack,
  Link,
  Typography,
  ListItemAvatar,
  Avatar,
} from '@mui/material';
import { Link as RouterLink } from 'react-router-dom';
import tittleToColor from 'shared/utils/stringColor';

// Generic card component to reuse the same structure for the
// different search card items as Institution, team, ...
function SearchCard({
  action, avatar, title, link,
  primary, secondary, tertiary,
  displayAvatar = false,
}) {
  return (
    <ListItem divider alignItems="flex-start" secondaryAction={action}>
      {displayAvatar
      && (
      <ListItemAvatar spacing={10}>
        <Avatar
          sx={{ bgcolor: tittleToColor(title), width: 45, height: 45 }}
          alt={title}
          src={avatar || 'dummy.png'}
        />
      </ListItemAvatar>
      )}
      <Stack direction="column" spacing={0.5}>
        {/* primary */}
        <Stack direction="row" spacing={2} alignItems="center">
          <Link
            sx={{ fontSize: '24px' }}
            component={RouterLink}
            to={link}
            underline="hover"
          >
            <b>{title}</b>
          </Link>
          {primary}
        </Stack>
        {/* secondary */}
        <Typography sx={{ color: 'text.secondary' }}>
          {secondary}
        </Typography>
        {/* tertiary */}
        <Stack
          direction="row"
          sx={{ color: 'text.secondary' }}
          alignItems="flex-end"
          spacing={2}
        >
          {tertiary}
        </Stack>
      </Stack>
    </ListItem>
  );
}

export default SearchCard;
