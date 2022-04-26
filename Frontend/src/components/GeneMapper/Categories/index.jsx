import {
  ExpandLess, ExpandMore, Visibility, VisibilityOff,
} from '@mui/icons-material';
import OpacityIcon from '@mui/icons-material/Opacity';
import {
  Box, Collapse, IconButton, ListItem, ListItemButton, Typography,
} from '@mui/material';
import React, { useState } from 'react';
import { colors } from 'shared/theme/colors';

const activatedColor = colors.primary['400'];
const deactivatedColor = colors.primary['200'];

function Category({
  title, values, colored, toggleColored,
}) {
  const [open, setOpen] = useState(false);

  return (
    <Box sx={{ display: 'flex', flexDirection: 'column' }}>
      <ListItemButton onClick={() => setOpen(!open)} sx={{ p: 0 }}>
        <ListItem
          secondaryAction={(
            <IconButton edge="end" onClick={(e) => { toggleColored(); e.stopPropagation(); }}>
              <OpacityIcon sx={{ color: colored ? activatedColor : deactivatedColor }} />
            </IconButton>
          )}
          sx={{ pl: 1 }}
        >
          {open
            ? <ExpandLess sx={{ transform: 'rotate(180deg)' }} />
            : <ExpandMore sx={{ transform: 'rotate(-90deg)' }} />}
          <Typography sx={{ pr: 1 }}>{title}</Typography>

        </ListItem>
      </ListItemButton>
      <Collapse in={open}>
        <Box sx={{ display: 'flex', flexDirection: 'column' }}>
          {values.map((value) => (
            <Value
              title={value.title}
              color={value.color}
              key={value.title}
            />
          ))}
        </Box>
      </Collapse>
    </Box>
  );
}

function Value({ title, color }) {
  const [visible, setVisible] = useState(true);

  return (
    <Box sx={{
      display: 'flex', alignItems: 'center', pl: 3, pr: 2,
    }}
    >
      <IconButton onClick={() => setVisible(!visible)}>
        { visible
          ? <Visibility sx={{ color }} />
          : <VisibilityOff sx={{ color: deactivatedColor }} />}
      </IconButton>
      <Typography sx={{ flexGrow: 1 }}>{title}</Typography>
    </Box>
  );
}

/**
 *
 * @param categories An object containing available categories and their values as key-value pairs
 */
function GeneMapperCategories({ categories }) {
  const [coloredCategoryTitle, setColoredCategoryTitle] = useState(null);

  return (
    <>
      {Object.keys(categories).map((title) => (
        <Category
          key={title}
          title={title}
          values={categories[title]}
          colored={title === coloredCategoryTitle}
          toggleColored={() => {
            setColoredCategoryTitle(title);
          }}
        />
      ))}
    </>

  );
}

export default GeneMapperCategories;
