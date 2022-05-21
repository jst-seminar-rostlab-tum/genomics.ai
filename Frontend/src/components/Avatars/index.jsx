import React from 'react';
import { AvatarGroup, Avatar } from '@mui/material';
import { grey } from '@mui/material/colors';
import MoreHorizIcon from '@mui/icons-material/MoreHoriz';
import stringToColor from 'shared/utils/stringColor';

// A component displaying a common group of avatrs
const Avatars = ({
  items,
  width = 23,
  height = 23,
  displayMoreIcon = true,
  // eslint-disable-next-line arrow-body-style
}) => {
  return (
    <AvatarGroup>
      {items.map((item) => (
        <Avatar
          key={item.src}
          sx={{ bgcolor: stringToColor(item.alt), width, height }}
          alt={item.alt}
          src={item.src || '/dummy.jpg'} // Not perfect but otherwise display default user profile icon as fallback rather then letter
        />
      ))}
      {displayMoreIcon && (
        <Avatar sx={{ bgcolor: grey[600], width, height }} src="/dummy.jpg">
          <MoreHorizIcon />
        </Avatar>
      )}
    </AvatarGroup>
  );
};

export default Avatars;
