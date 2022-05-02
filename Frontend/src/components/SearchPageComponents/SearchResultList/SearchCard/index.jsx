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

// Temporary function to generate different colors,
// change later to titleToColor(...) from feature/institutinOverview branch
function randomColor() {
  const hex = Math.floor(Math.random() * 0xffffff);
  return `#${hex.toString(16)}`;
}

// Generic card component to reuse the same structure for the
// different search card items as Institution, team, ...
function SearchCard({
  action, avatar, title, link,
  primary, secondary, tertiary,
}) {
  return (
    <ListItem divider alignItems="flex-start" secondaryAction={action}>
      <ListItemAvatar spacing={10}>
        <Avatar
          sx={{ bgcolor: randomColor(), width: 45, height: 45 }}
          alt={title}
          src={avatar || 'dummy.png'} // not nice, but fallback doesn't work properly and display default profile icon rather than the first letter of alt
        />
      </ListItemAvatar>
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
