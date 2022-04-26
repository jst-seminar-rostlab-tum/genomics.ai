import {
  Circle, ExpandLess, ExpandMore, FormatColorFill, FormatColorReset, Visibility, VisibilityOff,
} from '@mui/icons-material';
import {
  Box, Collapse, IconButton, ListItem, ListItemButton, Typography,
} from '@mui/material';
import React, { useState } from 'react';

function Category({
  title, values, colored, toggleColored,
}) {
  const [open, setOpen] = useState(false);

  return (
    <Box sx={{ display: 'flex', flexDirection: 'column' }}>
      <ListItemButton onClick={() => setOpen(!open)} sx={{ p: 0 }}>
        <ListItem
          secondaryAction={(
            <IconButton onClick={(e) => { toggleColored(); e.stopPropagation(); }}>
              {colored ? <FormatColorFill /> : <FormatColorReset />}
            </IconButton>
          )}
        >
          <Typography sx={{ pr: 1 }}>{title}</Typography>
          {open ? <ExpandLess /> : <ExpandMore />}

        </ListItem>
      </ListItemButton>
      <Collapse in={open}>
        <Box sx={{ display: 'flex', flexDirection: 'column' }}>
          {values.map((valueTitle) => <Value title={valueTitle} key={valueTitle} />)}
        </Box>
      </Collapse>
    </Box>
  );
}

function Value({ title }) {
  const [visible, setVisible] = useState(true);

  return (
    <Box sx={{ display: 'flex', alignItems: 'center', px: 3 }}>
      <IconButton edge="start" onClick={() => setVisible(!visible)}>
        {visible ? <Visibility /> : <VisibilityOff />}
      </IconButton>
      <Typography sx={{ flexGrow: 1, pr: 3 }}>{title}</Typography>
      <Circle />
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
